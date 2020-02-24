
    
rule call_compartment:
    input:
      cool = join(OUT, 'cool', '{sample}_{libtype}_' + f'{COMP_RES_STR}.cool'),
      gene_cov = join(TMP, 'gene_coverage.tsv')
    output: join(OUT, 'compartments', 'compartments_{sample}_{libtype}.bedgraph')
    params:
        comp_res = COMP_RES
    conda:
        'envs/hic_processing.yaml'
    script:
        "../scripts/compartment_utils.py"

rule merge_compartments:
  input : expand(join(OUT, 'compartments', 'compartments_{sample}_HiC.bedgraph'),sample=samples['name'])
  output: join(OUT, 'compartments', 'merged_compartments.tsv')
  run:
    # Run bin coordinates from ref file
    bedgraph = pd.read_csv(input[0], sep='\t', usecols=[0, 1, 2], names=['chrom', 'start', 'end'], header=None)
    # Append eigenvector of libraries iteratively
    for i, sample in enumerate(samples['name']):
        print(sample)
        bedgraph[sample] = pd.read_csv(input[i], sep='\t', usecols=[3])
    bedgraph.to_csv(output[0], index=False, sep='\t')


rule plot_eigenvectors:
    input: join(OUT, 'compartments', 'merged_compartments.tsv')
    output: join(OUT, 'figures', 'compartments', 'eigens', '{chrom}_eigen.svg')
    params:
      samples_df = samples.loc[:,:],
      chrom = lambda w: f"{w.chrom}"
    run:
      import matplotlib
      matplotlib.use('Agg')
      comp_df = pd.read_csv(input[0], sep='\t')
      plot_eigens(comp_df, params['samples_df'], chrom=params['chrom'], out=output[0])

# # Measure compartment similarity on whole genome using PCA
# rule plot_compartment_pca:
#   input: join(OUT, 'compartments', 'merged_compartments.tsv')
#   output: join(OUT, 'figures', 'compartments', 'compartments_pca_{variable}.pdf')
#   params:
#     samples_df = samples.loc[np.isin(samples.library, hic_libs), :],
#   run:
#     import matplotlib
#     matplotlib.use('Agg')
#     comp_df = pd.read_csv(input[0], sep='\t')
#     pca_compartments(comp_df, params['samples_df'], col_var=wildcards['variable'], out=output[0])