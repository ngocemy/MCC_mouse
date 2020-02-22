import os


def gene_coverage(refseq_file, chrom_size_file, bin_size=100000):

    # To return a bedgraph file with the gene coverage (genes per bin) for uscs downloaded annotation file.
    # refseqfile: path to USCS refererence file with annotations; chrom_size_file: chromosome size file also downloaded
    # from UCSC that will be chopped into bins

    # Separate genome into bins. In this example, we generate genomic bins at 100kb by default resolution by using bedtools:
    # Make a directory for gene density track outputs, for example:

    # Sort gene positions in order: descending or ascendingg?

    # Count the number of genes per genomic bin, plot over the reference genome and then export as a bigwig file:
    
    script = 
    """
    awk -vOFS='\t' '{ print $3,$5,$6}' ${refseq_file} > refseq_bed.txt
    bedtools makewindows -g ${chrom_size_file} -w ${bin_size} -s ${bin_size} > mm10.100kb.bin
    bedtools coverage -a mm10.100kb.bin -b mm10_refGene_sort.txt > data/tmp/gene_density.txt
    awk -vOFS='\t' '{ print $1,$2,$3,$7}' data/tmp/gene_density.txt > data/tmp/gene_coverage.tsv
    """
    os.system("bash -c '%s'" % script)

