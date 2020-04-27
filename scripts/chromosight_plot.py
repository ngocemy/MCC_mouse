import re
import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt

# Load detected patterns' tables
loops = pd.read_csv(snakemake.input["loops"], sep="\t")
borders = pd.read_csv(snakemake.input["borders"], sep="\t")

# Load Hi-C data in cool format
c = cooler.Cooler(snakemake.input["cool_file"])

# res = 2000

# # Load images (vignettes) around RAD21 interactions coordinates
# images_g = np.load(snakemake.input[0])
# images_m = np.load("quantify/rad21_metaphase/loops_quant.npy")

# # Load lists of RAD21 interactions coordinates with their loop scores
# # Compute loop size (i.e. anchor distance) for each RAD21 combination
# get_sizes = lambda df: np.abs(df.start2 - df.start1)
# loops_g = pd.read_csv("quantify/rad21_g1/loops_quant.txt", sep="\t")
# loops_g["loop_size"] = get_sizes(loops_g)
# loops_m = pd.read_csv("quantify/rad21_metaphase/loops_quant.txt", sep="\t")
# loops_m["loop_size"] = get_sizes(loops_m)

# # Merge data from both conditions into a single table
# loops_g["condition"] = "g1"
# loops_m["condition"] = "metaphase"
# loops_df = pd.concat([loops_g, loops_m]).reset_index(drop=True)
# images = np.concatenate([images_g, images_m])

# # Remove NaN scores (e.g. in repeated regions or overlap the matrix edge)
# nan_mask = ~np.isnan(loops_df["score"])
# loops_df = loops_df.loc[nan_mask, :]
# images = images[nan_mask, :, :]

# # The loop kernel can be loaded using chromosight.kernels.loops
# kernel = np.array(ck.loops["kernels"][0])
# pileup_kw = {"vmin": -1, "vmax": 1, "cmap": "seismic"}


# Select a region of interest
region = "chr5:55000000-56000000"
mat = c.matrix(sparse=False, balance=True).fetch(region)


def subset_region(df, region):
    """
    Given a pattern dataframe and UCSC region string, retrieve only patterns in that region.
    """
    # Split the region string at each occurence of - or : (yields 3 elements)
    chrom, start, end = re.split("[-:]", region)
    start, end = int(start), int(end)
    # Only keep patterns on the same chromosome as the region and
    # within the start-end interval
    subset = df.loc[
        (df.chrom1 == chrom)
        & (df.chrom2 == chrom)
        & (df.start1 >= start)
        & (df.start2 >= start)
        & (df.end1 < end)
        & (df.end2 < end),
        :,
    ]
    return subset


loops_sub = subset_region(loops, region)
borders_sub = subset_region(borders, region)
# hairpins_sub = subset_region(hairpins, region)

# Make genome-based bin numbers relative to the region
for df in [loops_sub, borders_sub]:
    df.bin1 -= c.extent(region)[0]
    df.bin2 -= c.extent(region)[0]

plt.imshow(mat ** 0.2, cmap="afmhot_r")
plt.scatter(
    loops_sub.bin2, loops_sub.bin1, edgecolors="blue", facecolors="none", label="loops"
)
plt.scatter(borders_sub.bin2, borders_sub.bin1, c="lightblue", label="borders")
# plt.scatter(hairpins_sub.bin2, hairpins_sub.bin1, c="green", label="hairpins")
plt.legend()
plt.savefig(snakemake.output[0])
