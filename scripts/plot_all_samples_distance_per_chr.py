import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import hicstuff.distance_law
from hicstuff.distance_law import slope_distance_law, import_distance_law, get_ylim
import numpy as np
import os.path

samples = snakemake.params["sample_name"]
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(12, 18), dpi=80)
fig.suptitle(
    f"Distance plot for MCC differentiation on chromosome:{snakemake.wildcards['chrom']}"
)
# Call the function plot_ps_slope(xs, ps, labels, fig_path=None, inf=3000, sup=None) from hicstuff.distancelaw
inf = 100000
sup = 10000000
xs, ps, labels = [], [], []
for (i, sample) in enumerate(snakemake.input[:]):
    sample_xs, sample_ps, sample_labels = import_distance_law(sample)
    xs.append(sample_xs[0])
    ps.append(sample_ps[0])
    labels.append(snakemake.params["sample_name"][i])
slope = slope_distance_law(xs, ps)

# Compute slopes from the curves
# Make the plot of distance law
# Give a range of color
cols = iter(cm.rainbow(np.linspace(0, 1, len(ps))))
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(18, 10))
plt.subplots_adjust(left=0.05, right=0.85, top=0.93, bottom=0.07)
ax1.set_xlabel("Distance (pb)", fontsize="x-large")
ax1.set_ylabel("P(s)", fontsize="x-large")
ax1.set_title("Distance law", fontsize="xx-large")
ylim = get_ylim(xs, ps, inf, sup)
ax1.set_ylim(0.9 * ylim[0], 1.1 * ylim[1])
for i in range(len(ps)):
    # Iterate on the different distance law array and take them by order of
    # size in order to have the color scale equivalent to the size scale
    col = next(cols)
    ax1.loglog(xs[i], ps[i], label=labels[i])
# Make the same plot with the slope
cols = iter(cm.rainbow(np.linspace(0, 1, len(slope))))
ax2.set_xlabel("Distance (pb)", fontsize="x-large")
ax2.set_ylabel("Slope", fontsize="x-large")
ax2.set_title("Slope of the distance law", fontsize="xx-large")
ax2.set_xlim([inf, sup])
ylim = get_ylim(xs, slope, inf, sup)
ax2.set_ylim(1.1 * ylim[0], 0.9 * ylim[1])
xs2 = [None] * len(xs)
for i in range(len(slope)):
    xs2[i] = xs[i][:-1]
    col = next(cols)
    ax2.semilogx(xs2[i], slope[i], label=labels[i], subsx=[2, 3, 4, 5, 6, 7, 8, 9])
ax2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.00), ncol=1, fontsize="large")


plt.savefig(snakemake.output[0])

