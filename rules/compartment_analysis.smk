
    
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

