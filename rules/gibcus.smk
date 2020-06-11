rule chromosight_detect:
    input: join(IN, 'gibcus2018','cool','{sample_gibcus}.cool')
    output: 
        txt = join(OUT,'gibcus2018', 'detect','{sample_gibcus}', '{pattern}_out.txt'),
        npy = join(OUT,'gibcus2018', 'detect','{sample_gibcus}', '{pattern}_out.npy')
    threads: 5
    shell:
        """
        chromosight detect -t {threads} --win-fmt npy --pattern={wildcards.pattern} {input} $(dirname {output.txt})
        """
rule chromosight_quantify:
    input: 
        txt = join(OUT,'gibcus2018', 'detect','{sample_gibcus}','{pattern}_out.txt'),
        cool = join(IN, 'gibcus2018','cool','{sample_gibcus}.cool')
    output: 
        txt = join(OUT,'gibcus2018', 'quantify','{sample_gibcus}', '{pattern}_quant.txt'),
        npy = join(OUT,'gibcus2018', 'quantify','{sample_gibcus}', '{pattern}_quant.npy')
    shell:
        """
        chromosight quantify --pattern={wildcards.pattern} \
                             --win-fmt npy \
                             {input.txt} \
                             {input.cool} \
                             $(dirname {output.txt})
        """
rule analyze_chromosight:
    input: 
        cool_file = join(IN, 'gibcus2018','cool','{sample_gibcus}.cool'),
        loops = join(OUT,'gibcus2018','detect','{sample_gibcus}','loops_out.txt'),
        borders = join(OUT,'gibcus2018','detect','{sample_gibcus}','borders_out.txt')
    output: join(OUT,'gibcus2018','plot', '{sample_gibcus}.svg')
    script:
        "../scripts/gibcus.py"
