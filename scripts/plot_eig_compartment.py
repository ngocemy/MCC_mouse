import cooltools
import cooltools.eigdecomp
import cooltools.expected
import cooltools.saddle
import bioframe
import matplotlib.pyplot as plt
import cooler
import numpy as np


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


def saddle_plot(eig, cool, out=None):

    c = cooler.Cooler(cool)
    bins = c.bins()[:]
    regions = [(chrom, 0, c.chromsizes[chrom]) for chrom in c.chromnames]
    # Remove NaN
    # eig = eig.dropna().reset_index(drop=True)
    # Digitize eigenvectors, i.e. group genomic bins into
    # equisized groups according to their eigenvector rank.
    Q_LO = 0.025  # ignore 2.5% of genomic bins with the lowest E1 values
    Q_HI = 0.975  # ignore 2.5% of genomic bins with the highest E1 values
    N_GROUPS = (
        38  # divide remaining 95% of the genome into 38 equisized groups, 2.5% each
    )
    q_edges = np.linspace(Q_LO, Q_HI, N_GROUPS + 1)
    # Filter track used for grouping genomic bins based on bins filtered out in Hi-C balancing weights
    # Doesn't do anything with eigenvectors from the same Hi-C data (hence commented out here),
    # but important for external data, such as ChIP-seq tracks
    # eig = cooltools.saddle.mask_bad_bins((cis_eigs[1], 'E1'), (c.bins()[:], 'weight'))

    # Calculate the lower and the upper values of E1 in each of 38 groups.
    group_E1_bounds = cooltools.saddle.quantile(eig["E1"], q_edges)
    # breakpoint()
    # Assign the group to each genomic bin according to its E1, i.e. "digitize" E1.
    digitized, hist = cooltools.saddle.digitize_track(
        group_E1_bounds, track=(eig, "E1"),
    )
    # Calculate the decay of contact frequency with distance (i.e. "expected")
    # for each chromosome.
    # Retrieve cooler file and extract bin table

    expected = cooltools.expected.cis_expected(c, regions, use_dask=True)

    # Make a function that returns observed/expected dense matrix of an arbitrary
    # region of the Hi-C map.
    get_matrix = cooltools.saddle.make_cis_obsexp_fetcher(c, (expected, "balanced.avg"))
    # Compute the saddle plot, i.e. the average observed/expected between genomic
    # ins as a function of their digitized E1.
    S, C = cooltools.saddle.make_saddle(
        get_matrix, group_E1_bounds, (digitized, "E1" + ".d"), contact_type="cis"
    )

    plt.imshow(
        np.log2(S / C)[1:-1, 1:-1], cmap="coolwarm", vmin=-1, vmax=1,
    )
    plt.colorbar(label="log2 obs/exp")
    if out is None:
        plt.show()
    else:
        plt.savefig(out)
