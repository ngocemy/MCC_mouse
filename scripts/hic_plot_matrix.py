# To create plots from cool files output (Hicstuff)
import matplotlib.pyplot as plt
import hicstuff.view as hcv
import cooler
from os.path import join
import numpy as np

df = snakemake.params["df"]
res = snakemake.params["resolution"]
reg = snakemake.params["region"]
mats = {}
for sample in df["name"]:
    c = cooler.Cooler(
        join(
            snakemake.params["out"],
            "mcool",
            sample + f"_normalized.mcool::/resolutions/{res}",
        )
    )
    mats[sample] = [
        c.matrix(sparse=False).fetch(snakemake.params["region"]),
        c.info["sum"],
    ]
    # sum is total number of contacts


# Plot multiple plots in one frame not using HiCstuff but imshow
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12, 18), dpi=80)
for i, (s, mat) in enumerate(mats.items()):
    a = ax.flatten()[i]
    a.set_aspect(1)  # make the plot symmetrical
    a.set_title(
        f"{s}_{df.loc[df['name'] == s, 'condition'].values[0]}_total: {mats[s][1]}"
    )  # Total no of contacts
    a.imshow(
        np.log(mats[s][0]), cmap="afmhot_r"
    )  # vmax np.percentile is to decrease the contrast
plt.suptitle(
    f"Normalized Hi-C matrix for MCC differentiation at region {reg} and resolution {res}"
)
plt.savefig(snakemake.output[0])
