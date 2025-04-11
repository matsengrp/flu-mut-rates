import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import permutations as perm
import os

sns.set_theme(font_scale=1.0, style="ticks", palette="colorblind")


all_count_paths = snakemake.input.all_count_paths
boxplot_paths = snakemake.output.boxplot_paths
mut_type_dirs = snakemake.params.mut_type_dirs


def make_plots(all_counts_path, boxplot_path, mut_type_dir):
    syn_counts_df = pd.read_csv(all_counts_path)
    syn_counts_df.query("is_synonymous == True", inplace=True)
    syn_counts_df.sort_values("mut_type", inplace=True)

    plot_counts_by_mut_type(syn_counts_df, boxplot_path)

    os.makedirs(mut_type_dir, exist_ok=True)
    for mut_type in map("".join, perm("AGCT", 2)):
        outpath = f"{mut_type_dir}{mut_type}.png"
        plot_counts_for_mut_type(syn_counts_df, mut_type, outpath)

    return None


def plot_counts_by_mut_type(counts_df, outpath):
    plt.figure(figsize=[8, 3])
    sns.boxplot(x="mut_type", y="counts", data=counts_df, whis=(5, 95))
    plt.yscale("log")
    sns.despine()
    plt.savefig(outpath)
    plt.close()
    return None


def plot_counts_for_mut_type(counts_df, mut_type, outpath):
    data = counts_df.query("mut_type==@mut_type")
    plt.figure(figsize=[12, 3])
    sns.scatterplot(x="site", y="counts", data=data, alpha=0.5)
    plt.title(f"{mut_type}, N={len(data)}")
    sns.despine()
    plt.savefig(outpath)
    plt.close()
    return None


any(make_plots(*p) for p in zip(all_count_paths, boxplot_paths, mut_type_dirs))
