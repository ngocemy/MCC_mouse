# MCC_mouse
This project is _in collaboration with Alice Meunier_ to investigate the changes in the higher order chromosome
architecture during the progression of MCC (multiciliated cells) differentiation by HiC analysis.
### Dependencies
The following python packages must be installed to run the analysis:
* python >= 3.7
* snakemake
* GNU awk
* anaconda
### Input file
MCC libraries are either from non-differentiated cells or differentiated cells which are stored as fastq files for the reverse (R2) and forward(R1)......*tobecontinued*
### Running the analysis
The analysis is written in python and use conda environments (or optionally, singularity containers) to manage dependencies. The analysis is centralized in a master `Snakefile` and the different parts of the analysis are written in modular workflows in the `rules` directory. Those workflows are called by the `Snakefile` when running the pipeline.

The snakefile reads a number of parameters from the file `config.yaml`, such as the path to the reference genome.
Samples metadata and path to fastq files are described in the files `samples.tsv` and `units.tsv`, respectively.


To run the snakefile on 8 CPUs, simply use:

```bash
snakemake -j 8 --use-conda
```

The workflow will run independently for each library, each run starting from the genome and a pairs file:

### Output files

Different 1D signals are computed from the pairs file:
* Read coverage
* Insulation score
* A/B compartments scores
These signals are then combined into a bedgraph file for each sample. Taking the example above, the output files should be named.
```
data
└──output
    ├── all_signals_lib1.bedgraph
    └── all_signals_lib2.bedgraph
```
*tobecontinue*
