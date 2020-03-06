# Display a region within a cool file and the associated patterns detected by chromosight
# cmdoret, 20200301

import pandas as pd
import cooler
import matplotlib.pyplot as plt

# Load input matrix and chromosight patterns
loops = pd.read_csv('demo_loops_axel/loops_small_out.txt', sep='\t')
c = cooler.Cooler('data/in/SRR6675327.cool')

# Define region of interest in UCSC format
chrom, start, end = 'chr19', 100000, 4000000
region = f'{chrom}:{start}-{end}'

# Get whole genome bins corresponding to region
start_bin, _ = c.extent(region)

# Load matrix
mat = c.matrix(balance=True, sparse=True).fetch(region)

# Subset loops in region
loops = loops.loc[
        (loops.chrom1 == chrom) & 
        (loops.start1 >= start) &
        (loops.end1 < end) &
        (loops.chrom2 == chrom) &
        (loops.start2 >= start) &
        (loops.end2 < end), :
]

# Make bins relative to zoom
loops.bin1 -= start_bin
loops.bin2 -= start_bin

# Plot everything
plt.imshow(mat.toarray() ** 0.2, cmap='afmhot_r')
plt.scatter(loops.bin2, loops.bin1, facecolors='none', edgecolors='blue')
plt.show()
