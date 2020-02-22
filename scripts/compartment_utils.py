from cooltools import eigdecomp


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
    track = pd.read_csv(track_file, sep="\t")
    # genecov = bioframe.tools.frac_gene_coverage(bins, "mm10")

    # Compute 3 eigen vectors and rank + correlate with gene density
    cis_vals, cis_eigs = eigdecomp.cooler_cis_eig(
        c,
        bins=track,
        n_eigs=3,
        phasing_track_col=track.columns[3],
        sort_metric="pearsonr",
    )

    bg = cis_eigs.loc[:, ["chrom", "start", "end", "E1", "E2", "E3"]]

    bg.to_csv(out_bedgraph, sep="\t", header=None, index=False, na_rep="nan")


compartment_to_bedgraph(
    snakemake.input["cool"], snakemake.input["gene_cov"], snakemake.output[0]
)

