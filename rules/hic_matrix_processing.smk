conda: "../envs/hic_processing.yaml"

rule matrix_from_hicstuff:
    input:
        fa=GENOME,
        fq1=lambda w: units.loc[units['name'] == w.sample, 'fq1'],
        fq2=lambda w: units.loc[units['name'] == w.sample, 'fq2']
    output: 
        out_dir = directory(join(TMP, 'hicstuff', '{sample}')),
        cool_file = join(OUT, 'cool', '{sample}','{sample}.cool')
    threads: 6
    log: "logs/matrix_from_hicstuff/{sample}.log"
    shell:
        """
        hicstuff pipeline -g {input.fa} \
                        -o {output.out_dir} \
                        -t {threads} \
                         {input.fq1} \
                        {input.fq2} \
                        -M 'cool'\
                        -n \
                        -e 10000 \
                        2> {log}
        mv {output.out_dir}/abs_fragments_contacts_weighted.cool {output.cool_file}

        """
rule multiresolution:
    input: join(OUT, 'cool', '{sample}','{sample}.cool')
    output: join(OUT, 'mcool','_'.join(['{sample}','normalized.mcool']))
    threads: 6
    shell: 
        """
        cooler zoomify \
                --balance \
                -o {output} \
                -r 10000,20000,50000,100000,200000 \
                -p {threads} \
                {input}
        """



rule plot_matrices:
    input: expand(join(OUT, 'mcool','{sample}_normalized.mcool'), sample=samples['name'])
    output: join(OUT,'plots', f'hic_plot_all_samples_{high_res}_{REGION}.svg')
    params:
        resolution=high_res,
        df = samples,
        out = OUT,
        region = REGION 
    script: "../scripts/hic_plot_matrix.py"