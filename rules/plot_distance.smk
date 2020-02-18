rule genenerating_distance_plot_matrice:
    input: join(OUT, 'cool','_'.join(['{sample}','normalized.mcool']))
    output: directory(join(OUT, 'plots','distance_plots',f'hic_plot_distance_from_matrix_{high_res}.svg))
    script: "../scripts/hic_distance_plot.py"

#Only keep the ++ or -- strand to avoid uncut and selfing

rule genenerating_distance_plot_pairs_file:
    input: 
        pairs = join(
            TMP, 'hicstuff','{sample}', 'tmp', '{sample}'+'.valid_idx.pairs'
        ),
        frags = join(
            TMP, 'hicstuff','{sample}', '{sample}'+'.frags.tsv'
        )
    output: join(OUT, 'distance_law', '{sample}'+'_plot_distance_from_pairs.txt')
    shell:
        """
        hicstuff distancelaw   --frags {input.frags} \
                     --inf 60000 \
                     --pairs {input.pairs}\
                     -O {output}

        """

rule split_distance_chr:
    input: join(OUT, 'distance_law', '{sample}'+'_plot_distance_from_pairs.txt')
    output: join(OUT, 'distance_law_split', '{sample}_{{chrom}}'+'_plot_distance_from_pairs.txt')
    shell: "grep {chrom} {input} > {output}"

rule plot_all_distance_per_chr:
    input: expand(join(OUT, 'distance_law_split', '{sample}_{{chrom}}'+'_plot_distance_from_pairs.txt'),chrom = chrID,sample=samples['name'])
    output: join(OUT, 'distance_law', 'plots','all_samples_{chrom}_plot_distance_from_pairs.pdf')
    params:
        sample_name = samples['name']
        
    script: 'plot_all_samples_distance_per_chr.py'

