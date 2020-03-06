rule chromosight_detect:
    input: join(IN, 'gibcus2018','cool','{filename}')
    output: directory(join(OUT,'gibcus2018', 'chromosight','{filename}'))
    shell:
        """
        chromosight detect -t6 --win-fmt npy {input} {output}
        chromosight detect --pattern=borders --win-fmt npy {input} {output}
        """
rule chromosight_quantify:
    input: 
        out_dir = join(OUT,'gibcus2018', 'chromosight','{filename}'),
        cool = join(IN, 'gibcus2018','cool','{filename}')
    output: directory(join(OUT,'gibcus2018', 'chromosight','quantify','{filename}'))
    shell:
        """
        chromosight quantify --pattern=loops {input[0]}/loops_out.txt {input[1]} {output}
        chromosight quantify --pattern=borders {input[0]}/borders_out.txt {input[1]} {output}
        """
rule analyze_chormosight:
    input:
    output:
    script:
        "../scripts/gibcus.py"
