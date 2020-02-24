from cooltools import eigdecomp
import cooler
import pandas as pd
import matplotlib.pyplot as plt


def compartment_to_bedgraph(cool, track_file, out_bedgraph):
    """
    Computes eigen vectors from a cooler file, correlate and ranks them according
    to an input track and writes a bedgraph file.
    Parameters
    ----------
    cool : str
        Path to the input cooler file. Can also be the URI in an mcool file. For
        example: ex.mcool::/resolutions/640000
    track_file : str
        Path to the bedgraph file with a track that should be positively correlated
        with compartment A. 4th column is used as the signal.
    bedgraph : str
        Path to the output bedgraph file with compartment info.
    """
    # Retrieve cooler file and extract bin table
    c = cooler.Cooler(cool)
    bins = c.bins()[:]
    track = pd.read_csv(
        track_file, sep="\t", names=["chrom", "start", "end", "percentage"]
    )
    # genecov = bioframe.tools.frac_gene_coverage(bins, "mm10")

    # Compute 3 eigen vectors and rank + correlate with gene density
    cis_vals, cis_eigs = eigdecomp.cooler_cis_eig(
        c,
        bins=track,
        n_eigs=3,
        phasing_track_col=track.columns[3],
        sort_metric="pearsonr",
    )

    bg = cis_eigs.loc[:, ["chrom", "start", "end", "E1"]]
    # cooler built-in function already calculate pearsonr correlation to find PC ~ AB and delete diagonal

    bg.to_csv(out_bedgraph, sep="\t", header=None, index=False, na_rep="nan")


# compartment_to_bedgraph(
#     snakemake.input["cool"], snakemake.input["gene_cov"], snakemake.output[0]
# )


def pca_compartment():
    """
    Do PCA analysis, return PC on PC1 of compartment analysis (PCA on PC)
    For analyzing inter-sample distance, or effect of condition on AB compartment
    --------------



    """


def plot_eigens(comp_df, samples_df, chrom=None, out=None):
    """
    Plot compartments of an input given chromosome for all samples
    
    Parameters
    ----------
    comp_df : pandas.DataFrame
        Table of eigenvectors. Each row represents a bin from the Hi-C matrix,
        columns should be 'chrom', 'start', 'end' followed by one column per
        sample, where headers match the sample names in samples_df.
    samples_df : pandas.DataFrame
        Table of samples. Each row represents a sample. Columns should at least
        include be 'library' (the sample name) and condition (2 values). 
    chrom : str
    out : bool
    """
    # Subset eigenvectors for the chromosome
    if chrom is not None:
        eigen = comp_df.loc[comp_df.chrom == chrom, :]
    else:
        eigen = comp_df
    # Remove NaN
    index_with_nan = eigen.index[eigen.isnull().any(axis=1)]
    eigen.drop(index_with_nan, 0, inplace=True)
    # Draw figure to contain 1 panel per sample (row)
    _, axes = plt.subplots(samples_df.shape[0], 1)
    # Remove top and bottom margins around panels
    plt.subplots_adjust(hspace=0)
    libs = samples_df.sort_values("condition").name
    pops = samples_df.sort_values("condition").condition
    # lcols = ["#1f77b4", "#ff7f0e"]
    # Draw panel for each sample
    for i, (lib, pop) in enumerate(zip(libs, pops)):
        # Draw two filled areas: one for each compartment
        compA = eigen[lib].copy()
        compA[compA <= 0] = 0
        compB = eigen[lib].copy()
        compB[compB > 0] = 0
        axes[i].fill_between(eigen["start"], 0, compA, color="g")
        axes[i].fill_between(eigen["start"], 0, compB, color="r")
        # Black horizontal axis at the center
        axes[i].axhline(0, c="black")
        # We don't need axis ticks eigenvectors are arbitrary, just show
        # the population of each sample
        axes[i].set_ylabel(f"{pop[0].upper()}\n{lib}")
        axes[i].set_yticks([])
        # Remove top frame line except for the first plot
        if i > 0:
            axes[i].spines["top"].set_visible(False)
        # Remove bottom frame line except for the last plot
        if i < (samples_df.shape[0] - 1):
            axes[i].set_xticks([])
            axes[i].spines["bottom"].set_visible(False)
        # Color ylabel according to condition
        # lcol = lcols[pop]
        # axes[i].yaxis.label.set_color(lcol)
    plt.suptitle(f"Compartments per individual chromosome, {chrom}")
    if out is None:
        plt.show()
    else:
        plt.savefig(out)
