rule filter:
    input:
        pairs = join(
            TMP, 'hicstuff','{sample}', 'tmp', ''.join(['{sample}','.valid_idx.pairs'])
        )
    output: 
        dir_out = directory(join(OUT,'filter','plots','{sample}')),
        pairs = join(TMP,'filter','plots','{sample}.pairs')
    shell:
        """
        hicstuff filter --plot \
                        -P {wildcards.sample} \
                        -f {output.dir_out} \
                        {input.pairs} \
                        {output.pairs}
        """