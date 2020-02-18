import matplotlib.pyplot as plt
import pandas as pd
import hicstuff.distancelaw as hhd
from hdd import plot_ps_slope, import_distance_law

# import pickle
samples = snakemake.params["sample_name"]


(fig,) = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(12, 18), dpi=80)
fig.suptitle("Distance plot for MCC differentiation on chromosome")
# Call the function plot_ps_slope(xs, ps, labels, fig_path=None, inf=3000, sup=None) from hicstuff.distancelaw

for i in range(1, 4):
    xs, ps, labels = import_distance_law(snakemake.input)
    ax[i,].set_aspect(1)  # make the plot symmetrical

    ax[i,].set_title(samples[i])
    ax[i,].plot_ps_slope(xs, ps, labels, inf=60000)
# for a in ax.flat:
#     a.set(xlabel='Genomic distance (bp)', ylabel='Number of events or contacts')

# # Hide x labels and tick labels for top plots and y ticks for right plots.
# for a in ax.flat:
#     a.label_outer()

plt.savefig(snakemake.output[0])

