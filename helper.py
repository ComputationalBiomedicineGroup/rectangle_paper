import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
import numpy as np

def match_dfs2(df1, df2):
    df1 = df1.loc[:, df2.columns]
    df1 = df1.loc[df1.index.isin(df2.index)]
    df2 = df2.loc[df2.index.isin(df1.index)]
    df1 = df1.sort_index()
    df2 =  df2.sort_index()
    return df1, df2

def match_dfs(df1, df2):
    common_columns = list(set(df1.columns).intersection(df2.columns))
    df1_common = df1[common_columns]
    df2_common = df2[common_columns]
    df1_common = df1_common.loc[df1_common.index.isin(df2_common.index)]
    df2_common = df2_common.loc[df2_common.index.isin(df1_common.index)]
    df1_common = df1_common.sort_index()
    df2_common =  df2_common.sort_index()
    return df1_common, df2_common







def create_evaluation_df(expected_df, estimated_df):
    evalutation_df = expected_df.T.melt()
    estimated_df_e = estimated_df.T.melt()
    evalutation_df['result'] = estimated_df_e['value']
    evalutation_df.columns = ['Cell type','True Fractions' , 'Estimated Fractions']
    return evalutation_df

def evaulate_total_df(df):
    rsme_total = np.sqrt(np.mean((df['Estimated Fractions'] - df['True Fractions']) ** 2))
    r_value_total, p_value_total =  pearsonr(df['Estimated Fractions'], df['True Fractions'])
    return rsme_total, r_value_total, p_value_total


def evaluate_by_cell_type(gold_std, estimated):
    rsmes = []
    r_values = []
    p_values = []

    for i in range(0, len(gold_std)):
        pearson_r = pearsonr(estimated.iloc[i, :], gold_std.iloc[i, :])
        r_values.append(pearson_r[0])
        p_values.append(pearson_r[1])
        rsmes.append(np.sqrt(np.mean((estimated.iloc[i, :] - gold_std.iloc[i, :]) ** 2)))
    return rsmes, p_values, r_values

def create_scatter_plot_total(df, rsme, p_value, r_value, title="Total", y_axis_label="Estimated Fractions", x_axis_label="True Fractions", save_path=None, hue="Cell type", size=None):
    cell_palette = {
        'Brca': '#66C2A5',
        'Luad': '#FC8D62',
        'Lusc': '#8DA0CB',
        'CPE': '#FC8D62',
        'B cells': '#66C2A5',
        'Macrophages': '#FC8D62',
        'mDC': '#8DA0CB',
        'Monocytes': '#E78AC3',
        'NK cells': '#DDD700',
        'T cells': '#FFD92F',
        'T cells CD4 conv': '#B3B3B3',
        'T cells CD4 Total': '#B3B3B3',
        'T cells CD8': '#A6D854',
        'Tregs': '#E5C494',
        'pDC': '#E5C494',
        'Neutrophils': '#B3B3B3',
        'Stromal cells': '#66C2A5',
        'Plasma cells': '#FC8D62',
        'ILC': '#8DA0CB',
        'Platelet': '#E78AC3'  # Assuming you want to keep this color as it was not provided in the new list
    }
    if size:
        plt.figure(figsize=size)
    # Create a scatter plot with different colored points
    g = sns.scatterplot(x="True Fractions", y="Estimated Fractions",hue=hue, data=df, palette=cell_palette)

    sns.regplot(x="True Fractions", y="Estimated Fractions", data=df, scatter=False, color='#444444', ci=None)
    text = f"r = {round(r_value, 3)}, p = {round(p_value, 10)} \nRMSE = {round(rsme, 3)}\nn = {len(df)}"



    # label the axes
    g.set(xlabel=x_axis_label, ylabel=y_axis_label)

    g.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    xlabel = g.set_xlabel(x_axis_label)
    ylabel = g.set_ylabel(y_axis_label)

    # Make the labels bigger and bold
    xlabel.set_size(12)
    xlabel.set_weight('bold')
    ylabel.set_size(12)
    ylabel.set_weight('bold')


    g.set_title(title, fontsize=13, weight='bold')
    g.annotate(text, xy=(0.01, 0.95), xycoords="axes fraction",
                xytext=(5, -5), textcoords="offset points",
                ha="left", va="top")

    # Set xlim and ylim to start from (0, 0)
    g.set_xlim(0, g.get_xlim()[1])
    g.set_ylim(0, g.get_ylim()[1])



    # Add diagonal dashed line
    xlim = g.get_xlim()[1]
    ylim = g.get_ylim()[1]


    diag_line, = g.plot([0, min(xlim, ylim)], [0, min(xlim, ylim)], ls='--', color='gray')
    if save_path:
        g.get_figure().savefig(save_path, format='pdf', bbox_inches='tight')
    return g


def create_scatter_plot_total_comparison(df1, rsme1, p_value1, r_value1, df2, rsme2, p_value2, r_value2, title1="Graph 1 Title", title2="Graph 2 Title", save_path=None):
    # Create two subplots side by side
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot the first scatter plot on the left
    ax1 = sns.scatterplot(x="True Fractions", y="Estimated Fractions", hue="Cell type", data=df1, ax=axes[0])
    sns.regplot(x="True Fractions", y="Estimated Fractions", data=df1, scatter=False, color='#444444', ci=None, ax=axes[0])
    text1 = f"r = {round(r_value1, 3)}, p = {round(p_value1, 10)}\nRMSE = {round(rsme1, 3)}\nn = {len(df1)}"
    ax1.set_title(title1)
    ax1.annotate(text1, xy=(0.01, 0.95), xycoords="axes fraction",
                 xytext=(5, -5), textcoords="offset points",
                 ha="left", va="top")
    ax1.set_xlim(0, ax1.get_xlim()[1])
    ax1.set_ylim(0, ax1.get_ylim()[1])
    xlim1 = ax1.get_xlim()[1]
    ylim1 = ax1.get_ylim()[1]
    diag_line1, = ax1.plot([0, min(xlim1, ylim1)], [0, min(xlim1, ylim1)], ls='--', color='gray')
    ax1.legend_.remove()
    # Plot the second scatter plot on the right
    ax2 = sns.scatterplot(x="True Fractions", y="Estimated Fractions", hue="Cell type", data=df2, ax=axes[1])
    sns.regplot(x="True Fractions", y="Estimated Fractions", data=df2, scatter=False, color='#444444', ci=None, ax=axes[1])
    text2 = f"r = {round(r_value2, 3)}, p = {round(p_value2, 10)}\nRMSE = {round(rsme2, 3)}\nn = {len(df2)}"
    ax2.set_title(title2)
    ax2.annotate(text2, xy=(0.01, 0.95), xycoords="axes fraction",
                 xytext=(5, -5), textcoords="offset points",
                 ha="left", va="top")
    ax2.set_xlim(0, ax2.get_xlim()[1])
    ax2.set_ylim(0, ax2.get_ylim()[1])
    ax2.legend_.remove()
    xlim2 = ax2.get_xlim()[1]
    ylim2 = ax2.get_ylim()[1]
    diag_line2, = ax2.plot([0, min(xlim2, ylim2)], [0, min(xlim2, ylim2)], ls='--', color='gray')

    # Add labels to the subplots
    axes[0].set_xlabel("True Fractions")
    axes[0].set_ylabel("Estimated Fractions")
    axes[1].set_xlabel("True Fractions")
    axes[1].set_ylabel("Estimated Fractions")

    # Create a single legend for both subplots
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5))

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure in PDF format if a save_path is provided
    if save_path:
        fig.savefig(save_path, format='pdf', bbox_inches='tight')





def create_scatter_plot_celltypes(df, rsmes, p_values, r_values, save_path=None):
    g = sns.lmplot(
    data=df, x="True Fractions", y="Estimated Fractions", col="Cell type", hue="Cell type", palette="colorblind",col_wrap=3, ci=None, scatter_kws={"s": 10, "alpha": 1}, facet_kws={"sharex": False, "sharey": False}
    )

    i = 0
    for row_val, ax in g.axes_dict.items():
        text = f"r = {round(r_values[i], 3)}, p = {round(p_values[i], 5)} \nRSME = {round(rsmes[i], 3)}\nn = {int(len(df)/ len(p_values))}"
        ax.annotate(text, xy=(0.01,0.95), xycoords="axes fraction",
                        xytext=(5,-5), textcoords="offset points",
                        ha="left", va="top")
        ax.set_xlim(0, ax.get_xlim()[1])
        ax.set_ylim(0, ax.get_ylim()[1])

        # Add diagonal dashed line
        xlim = ax.get_xlim()[1]
        ylim = ax.get_ylim()[1]
        diag_line, = ax.plot([0, min(xlim, ylim)], [0, min(xlim, ylim)], ls='--', color='gray')
        i = i + 1
    if save_path:
        g.savefig(save_path, format='pdf', bbox_inches='tight')
    g
