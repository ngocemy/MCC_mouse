conda: "../envs/hic_processing.yaml"

rule chrom_sizes:
  input: GENOME
  output: join(TMP, "chrom.sizes")
  params:
    genome = GENOME,
    res = MAX_RES,
    digest_dir = join(TMP, 'digest')
  singularity: 'docker://koszullab/hicstuff'

  shell:
    """
    hicstuff digest -e {params.res} \
                    -o {params.digest_dir} \
                    {input}
    tail -n +2 {params.digest_dir}/info_contigs.txt \
      | awk -vOFS='\t' '{{print $1,$2}}' \
      > {output}
    """

# rule matrix_from_hicstuff:
#     input:
#         fa=GENOME,
#         fq1=lambda w: units.loc[units['name'] == w.sample, 'fq1'],
#         fq2=lambda w: units.loc[units['name'] == w.sample, 'fq2']
#     output: 
#         out_dir = directory(join(TMP, 'hicstuff', '{sample}')),
#         cool_file = join(OUT, 'cool', '{sample}','{sample}.cool')
#     threads: 6
#     log: "logs/matrix_from_hicstuff/{sample}.log"
#     shell:
#         """
#         hicstuff pipeline -g {input.fa} \
#                         -o {output.out_dir} \
#                         -t {threads} \
#                          {input.fq1} \
#                         {input.fq2} \
#                         -M 'cool'\
#                         -n \
#                         -e 10000 \
#                         2> {log}
#         mv {output.out_dir}/abs_fragments_contacts_weighted.cool {output.cool_file}

#         """
rule matrix_compartment: # using pairs files to convert to cool for compartment analysis:
  input: 
    pairs = join(TMP, 'hicstuff','{sample}', 'tmp', 'valid.pairs'),
    chromsizes = join(TMP, 'chrom.sizes')
  output:
    cool_comp = join(OUT, 'cool', '{sample}_{libtype}_' + f'{COMP_RES_STR}.cool')
  params:
    comp_res = COMP_RES
  threads: 6
  shell: 
    """
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
                       {input.chromsizes}:{params.comp_res} \
                       {input.pairs} {output.cool_comp}
    cooler balance -p {threads} --mad-max 5 {output.cool_comp}
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
  output: join(OUT,'plots', f'hic_plot_all_samples_{LOW_RES}_{REGION}.svg')
  params:
    resolution=LOW_RES,
    df = samples,
    out = OUT,
    region = REGION 
  script: "../scripts/hic_plot_matrix.py"


rule plot_matrices_all_chr_each_sample:
  input: join(OUT, 'cool', '{sample}_HiC_' + f'{COMP_RES_STR}.cool')
  output: join(OUT, 'plots', 'all_chr','hic_plot_all_chrom_{sample}_' +  f'{LOW_RES_STR}.svg')
  params:
    resolution = LOW_RES,
    sample = lambda w: f"{w.sample}",
  run:
    res = params['resolution']
    res = bp_to_suffix(res)
    c = cooler.Cooler(input[0])
    sample = params['sample']
    plot_matrix(c, chrID,sample,res,output[0])

def filter_pcr_dup(pairs_idx_file, filtered_file):
    """
    Filter out PCR duplicates from a coordinate-sorted pairs file using
    overrrepresented exact coordinates. If multiple fragments have two reads
    with the exact same coordinates, only one of those fragments is kept.
    Parameters
    ----------
    pairs_idx_file : str
        Path to an indexed pairs file containing the Hi-C reads.
    filtered_file : str
        Path to the output pairs file after removing duplicates.
    """
    # Keep count of how many reads are filtered
    filter_count = 0
    reads_count = 0
    # Store header lines
    header = hio.get_pairs_header(pairs_idx_file)
    with open(pairs_idx_file, "r") as pairs, open(filtered_file, "w") as filtered:
        # Copy header lines to filtered file
        for head_line in header:
            filtered.write(head_line + "\n")
            next(pairs)

        # Use csv methods to easily access columns
        paircols = [
            "readID",
            "chr1",
            "pos1",
            "chr2",
            "pos2",
            "strand1",
            "strand2",
            "frag1",
            "frag2",
        ]
        # Columns used for comparison of coordinates
        coord_cols = [col for col in paircols if col != "readID"]
        pair_reader = csv.DictReader(pairs, delimiter="\t", fieldnames=paircols)
        filt_writer = csv.DictWriter(filtered, delimiter="\t", fieldnames=paircols)

        # Initialise a variable to store coordinates of reads in previous pair
        prev = {k: 0 for k in paircols}
        for pair in pair_reader:
            reads_count += 1
            # If coordinates are the same as before, skip pair
            if all(pair[pair_var] == prev[pair_var] for pair_var in coord_cols):
                filter_count += 1
                continue
            # Else write pair and store new coordinates as previous
            else:
                filt_writer.writerow(pair)
                prev = pair
        logger.info(
            "%d%% PCR duplicates have been filtered out (%d / %d pairs) "
            % (100 * round(filter_count / reads_count, 3), filter_count, reads_count)
        )