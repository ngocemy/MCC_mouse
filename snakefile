#!/bin/env snakemake -s
# This file can be run using snakemake. It runs all different HBV cancer 
# analyses on the input pairs files
# cmdoret, 20190501

from snakemake.utils import validate
from os.path import join
import numpy as np
import pandas as pd
import sys
from glob import glob

### USER DEFINED VARIABLES
# =============================================================================

# Define samples on which the code should run
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t', dtype=str, comment='#')
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep='\t', dtype=str, comment='#')
# Enforce str in index
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

# Set input / output paths
DATA_DIR = 'data'
IN = join(DATA_DIR, 'input')
OUT = join(DATA_DIR, 'output')
TMP = join(DATA_DIR, 'tmp')
GENOME = config['reference']
high_res = config['high_res']
low_res = config['low_res']
REGION = config['region']
chrID = ["chr{}".format(x) for x in list(range(1, 19)) + ["X"]]

wildcard_constraints:
  sample="|".join(samples['name'])


#breakpoint()
include : 'rules/hic_matrix_processing.smk'
#include : 'rules/plot_distance.smk'
#include : 'rules/filter_uncut_loop.smk'

rule all:
  input: 
    expand(join(TMP, 'hicstuff', '{sample}'),sample=samples['name']),
    #expand(join(TMP, 'hicstuff', '{sample}','.'.join(['{sample}','cool'])),sample=samples['name'])
    join(OUT,'plots', f'hic_plot_all_samples_{high_res}_{REGION}.svg')
    #expand(join(OUT,'filter','plots','{sample}'),sample=samples['name'])
    #expand(join(OUT, 'distance_law', 'plots','all_samples_{chrom}_plot_distance_from_pairs.pdf'),chrom=chrID)
