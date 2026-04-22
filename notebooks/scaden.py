#!/usr/bin/env python3

from __future__ import annotations

import gc
import os
import platform
import subprocess
import threading
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

MIB = 1024 * 1024
GIB = 1024 * 1024 * 1024


def _parse_mem_to_bytes(mem_str: str) -> int:
    # docker stats format: "123.4MiB / 7.7GiB"
    raw = mem_str.split("/")[0].strip()
    if not raw:
        return 0
    num = ""
    unit = ""
    for ch in raw:
        if ch.isdigit() or ch == ".":
            num += ch
        elif not ch.isspace():
            unit += ch
    if not num:
        return 0
    value = float(num)
    unit = unit.strip().lower()
    factors = {
        "b": 1,
        "kb": 1000,
        "mb": 1000**2,
        "gb": 1000**3,
        "kib": 1024,
        "mib": 1024**2,
        "gib": 1024**3,
    }
    return int(value * factors.get(unit, 1))


class DockerPeakMemory:
    def __init__(self, container_name: str, interval: float = 0.5) -> None:
        self.container_name = container_name
        self.interval = interval
        self.max_bytes = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def start(self) -> None:
        self.max_bytes = 0
        self._stop.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self) -> None:
        while not self._stop.is_set():
            res = subprocess.run(
                ["docker", "stats", "--no-stream", "--format", "{{.MemUsage}}", self.container_name],
                capture_output=True,
                text=True,
            )
            if res.returncode == 0 and res.stdout.strip():
                self.max_bytes = max(self.max_bytes, _parse_mem_to_bytes(res.stdout.strip()))
            time.sleep(self.interval)

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join()


def _require_tokens(template: str, tokens: list[str]) -> None:
    missing = [tok for tok in tokens if tok not in template]
    if missing:
        raise ValueError(
            "SCADEN_CMD_TEMPLATE is missing placeholders: "
            + ", ".join(missing)
            + ". Required: {MIXTURE}, {SC_H5AD_DIR}, {OUTDIR}."
        )


def _stratified_obs_sample(adata: anndata.AnnData, max_cells: int, seed: int) -> anndata.AnnData:
    if adata.n_obs <= max_cells:
        return adata

    rng = np.random.default_rng(seed)
    labels = adata.obs["cell_type"].astype(str)
    label_counts = labels.value_counts()
    total = int(labels.shape[0])
    picked: list[str] = []

    for label, count in label_counts.items():
        n_pick = max(1, int(round(max_cells * (count / total))))
        obs_names = labels.index[labels == label].to_numpy()
        n_pick = min(n_pick, obs_names.shape[0])
        picked.extend(rng.choice(obs_names, size=n_pick, replace=False).tolist())

    if len(picked) > max_cells:
        picked = rng.choice(np.array(picked), size=max_cells, replace=False).tolist()
    elif len(picked) < max_cells:
        remaining = np.setdiff1d(adata.obs_names.to_numpy(), np.array(picked))
        need = min(max_cells - len(picked), remaining.shape[0])
        if need > 0:
            picked.extend(rng.choice(remaining, size=need, replace=False).tolist())

    return adata[picked, :].copy()


def main() -> None:
    notebooks_dir = Path(__file__).resolve().parent
    repo_root = notebooks_dir.parent
    data_dir = repo_root / "data"
    fig3_dir = data_dir / "fig3"
    fig3_dir.mkdir(parents=True, exist_ok=True)

    adata_path = data_dir / "hao_adata.h5ad"
    bulks_path = data_dir / "bulks_finotello_tpm.csv"
    if not adata_path.exists():
        raise FileNotFoundError(f"Missing input file: {adata_path}")
    if not bulks_path.exists():
        raise FileNotFoundError(f"Missing input file: {bulks_path}")

    adata_hao = anndata.read_h5ad(adata_path)
    if "cell_type" not in adata_hao.obs:
        raise KeyError("`cell_type` column was not found in hao_adata.h5ad .obs")

    bulks_fino = pd.read_csv(bulks_path, index_col=0)
    sc_counts = adata_hao.to_df().T  # genes x cells

    common_genes = bulks_fino.index.intersection(sc_counts.index)
    if len(common_genes) == 0:
        raise ValueError("No overlapping genes between single-cell data and bulk data.")
    bulks_common = bulks_fino.loc[common_genes]

    mixture_path = fig3_dir / "scaden_mixture_full.tsv"
    bulks_common.to_csv(mixture_path, sep="\t")

    sim_max_cells = 90000
    print(f"Simulation input cap: {sim_max_cells} cells")
    sim_seed = int(os.environ.get("SCADEN_SIM_SEED", "42"))
    adata_for_sim = adata_hao if sim_max_cells <= 0 else _stratified_obs_sample(adata_hao, sim_max_cells, sim_seed)
    suffix = "full" if sim_max_cells <= 0 else f"max{adata_for_sim.n_obs}"
    sim_n_samples = 1000
    sim_cells_per_bulk = 100
    train_steps = 5000

    scaden_input_dir = fig3_dir / "scaden_h5ad_input" / suffix
    scaden_input_dir.mkdir(parents=True, exist_ok=True)
    sc_h5ad_path = scaden_input_dir / "hao_counts.h5ad"
    if not sc_h5ad_path.exists():
        sc_adata = adata_for_sim.copy()
        sc_adata.obs = sc_adata.obs.copy()
        sc_adata.obs["Celltype"] = sc_adata.obs["cell_type"].astype(str).values
        if not sparse.issparse(sc_adata.X):
            sc_adata.X = sparse.csr_matrix(sc_adata.X)
        sc_adata.write_h5ad(sc_h5ad_path)

    inner_cmd_template = os.environ.get("SCADEN_CMD_TEMPLATE", "").strip()
    if not inner_cmd_template:
        inner_cmd_template = (
            "scaden simulate --data {SC_H5AD_DIR} --data-format h5ad --pattern '*_counts.h5ad' "
            f"--out {{OUTDIR}}/simulate -n {sim_n_samples} --cells {sim_cells_per_bulk} --prefix {{OUTDIR}}/simulate/data && "
            "scaden process {OUTDIR}/simulate/data.h5ad {MIXTURE} --processed_path {OUTDIR}/processed.h5ad && "
            f"scaden train {{OUTDIR}}/processed.h5ad --model_dir {{OUTDIR}}/model --steps {train_steps} && "
            "scaden predict {MIXTURE} --model_dir {OUTDIR}/model --outname {OUTDIR}/predictions.txt"
        )
    _require_tokens(inner_cmd_template, ["{MIXTURE}", "{SC_H5AD_DIR}", "{OUTDIR}"])
    scaden_image = os.environ.get("SCADEN_DOCKER_IMAGE", "rectangle-scaden:arm64").strip()
    if not scaden_image:
        raise ValueError("SCADEN_DOCKER_IMAGE is empty.")
    docker_platform = os.environ.get("SCADEN_DOCKER_PLATFORM", "linux/arm64/v8").strip()
    docker_stats_interval = float(os.environ.get("SCADEN_DOCKER_STATS_INTERVAL", "0.5"))

    n_runs = 1
    runtime_records: list[dict[str, float | int]] = []
    memory_records: list[dict[str, float | int]] = []

    n_cells = int(adata_hao.n_obs)
    print(f"Benchmarking Scaden on full Hao dataset with {n_cells} cells")
    print(f"Simulation/training input cells: {adata_for_sim.n_obs}")
    print(
        "Scaden simulate/train params: "
        f"n_samples={sim_n_samples}, cells_per_bulk={sim_cells_per_bulk}, steps={train_steps}"
    )

    for run in range(1, n_runs + 1):
        run_outdir = fig3_dir / f"scaden_full_run_{run}"
        run_outdir.mkdir(parents=True, exist_ok=True)
        (run_outdir / "simulate").mkdir(parents=True, exist_ok=True)

        mixture_in_container = f"/workspace/{mixture_path.relative_to(repo_root).as_posix()}"
        sc_h5ad_in_container = f"/workspace/{sc_h5ad_path.relative_to(repo_root).as_posix()}"
        sc_h5ad_dir_in_container = f"/workspace/{scaden_input_dir.relative_to(repo_root).as_posix()}"
        outdir_in_container = f"/workspace/{run_outdir.relative_to(repo_root).as_posix()}"

        inner_cmd = (
            inner_cmd_template.replace("{MIXTURE}", mixture_in_container)
            .replace("{SC_H5AD}", sc_h5ad_in_container)
            .replace("{SC_H5AD_DIR}", sc_h5ad_dir_in_container)
            .replace("{OUTDIR}", outdir_in_container)
        )
        container_name = f"scaden-bench-{run}-{int(time.time())}"
        docker_cmd = ["docker", "run", "--rm", "--name", container_name]
        if docker_platform:
            docker_cmd.extend(["--platform", docker_platform])
        docker_cmd.extend(
            [
                "-v",
                f"{repo_root}:/workspace",
                "-w",
                "/workspace",
                scaden_image,
                "/bin/bash",
                "-lc",
                inner_cmd,
            ]
        )

        print(f"Run {run}/{n_runs}")
        gc.collect()
        time.sleep(0.2)

        sampler = DockerPeakMemory(container_name=container_name, interval=docker_stats_interval)
        sampler.start()
        t0 = time.perf_counter()
        try:
            subprocess.run(docker_cmd, check=True)
        except subprocess.CalledProcessError as exc:
            if exc.returncode == 132:
                host_arch = platform.machine().lower()
                raise RuntimeError(
                    "Scaden container crashed with exit code 132 (illegal instruction), "
                    "which usually means CPU architecture incompatibility (e.g. amd64 image on arm64 host). "
                    f"Current host arch: {host_arch}. "
                    "Use an ARM-compatible Scaden image via SCADEN_DOCKER_IMAGE "
                    "(default in this repo is rectangle-scaden:arm64), "
                    "or explicitly set SCADEN_DOCKER_PLATFORM."
                ) from exc
            if exc.returncode == 137:
                raise RuntimeError(
                    "Scaden container was OOM-killed (exit 137), typically during `scaden simulate` loading "
                ) from exc
            raise
        finally:
            t1 = time.perf_counter()
            sampler.stop()

        elapsed_s = t1 - t0
        peak_mem_mib = sampler.max_bytes / MIB
        peak_mem_gib = sampler.max_bytes / GIB

        runtime_records.append({"n_cells": n_cells, "run": run, "time": elapsed_s})
        memory_records.append({"n_cells": n_cells, "run": run, "memory_mib": peak_mem_mib, "memory_gib": peak_mem_gib})
        print(f"  {elapsed_s:.4f}s, {peak_mem_mib:.2f} MiB ({peak_mem_gib:.3f} GiB)")

    runtime_df = pd.DataFrame(runtime_records)
    memory_df = pd.DataFrame(memory_records)

    runtime_file = fig3_dir / "scaden_runtime_full_py.csv"
    memory_file = fig3_dir / "scaden_memory_full_py.csv"
    runtime_df.to_csv(runtime_file, index=False)
    memory_df.to_csv(memory_file, index=False)

    runtime_summary = runtime_df["time"].agg(["mean", "std"])
    memory_summary_mib = memory_df["memory_mib"].agg(["mean", "std"])
    memory_summary_gib = memory_df["memory_gib"].agg(["mean", "std"])
    print("\n--- Scaden Full-Dataset Benchmark Summary ---")
    print(f"Runtime mean={runtime_summary['mean']:.4f}s sd={runtime_summary['std']:.4f}")
    print(f"Memory mean={memory_summary_mib['mean']:.2f} MiB sd={memory_summary_mib['std']:.2f}")
    print(f"Memory mean={memory_summary_gib['mean']:.3f} GiB sd={memory_summary_gib['std']:.3f}")
    print(f"Saved runtime: {runtime_file}")
    print(f"Saved memory: {memory_file}")


if __name__ == "__main__":
    main()
