#!/bin/env snakemake -s
# This file can be run using snakemake. It runs all different HBV cancer 
# analyses on the input pairs files


from snakemake.utils import validate
from os.path import join
import numpy as np
import pandas as pd
import sys
from glob import glob

def bp_to_suffix(size):
    """
    Given a number of basepairs, returns the notation with suffix
    
    Examples
    --------
    >>> get_bp_scale(10000)
    "10kb"
    """
    # Defind mapping between powers of 10 and suffixes
    pow_to_suffix = {0: "bp", 3: "kb", 6: "Mb", 9: "Gb", 12: "Tb"}
    sorted_pows = sorted(pow_to_suffix.keys())
    # Find order of magnitude of input size
    input_len_pow = int(np.log10(size))
    # Find which power matches order of magnitude
    valid_pow_idx = max(0, np.searchsorted(sorted_pows, input_len_pow, side='right') - 1)
    input_valid_pow = sorted_pows[valid_pow_idx]
    # Get corresponding suffix
    suffix = pow_to_suffix[input_valid_pow]
    scale = 10 ** input_valid_pow
    str_bp = f"{int(size // scale)}{suffix}"
    return str_bp



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
COMP_RES = config['contact_maps']['comp_res']
COMP_RES_STR = bp_to_suffix(COMP_RES)
MAX_RES = config['contact_maps']['max_res']
wildcard_constraints:
  sample="|".join(samples['name'])
  #chrom="|".join(chrID)

#breakpoint()
#include : 'rules/hic_matrix_processing.smk'
#include : 'rules/plot_distance.smk'
#include : 'rules/filter_uncut_loop.smk'
include: 'rules/compartment_analysis.smk'
rule all:
  input:  
    #expand(join(OUT, 'cool', '{sample}_{libtype}_' + f'{COMP_RES_STR}.cool'),sample=samples['name'],libtype=units['lib_type'])  
    #expand(join(OUT, 'distance_law', '{sample}'+'_plot_distance_from_pairs.txt'),sample=samples['name'])
    #expand(join(TMP, 'hicstuff', '{sample}'),sample=samples['name']),
    #expand(join(TMP, 'hicstuff', '{sample}','.'.join(['{sample}','cool'])),sample=samples['name'])
    #join(OUT,'plots', f'hic_plot_all_samples_{high_res}_{REGION}.svg')
    #expand(join(OUT,'filter','plots','{sample}'),sample=samples['name'])
    #expand(join(OUT, 'distance_law', 'plots','all_samples_{chrom}_plot_distance_from_pairs.pdf'),chrom=chrID)
    expand(join(OUT, 'compartments', 'compartments_{sample}_{libtype}.bedgraph'),sample=samples['name'],libtype=units['lib_type'])
