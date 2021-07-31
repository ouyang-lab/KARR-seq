# KARR-seq

KARR-seq is a method that reveals RNA-RNA interactions in a transcriptome-wide fashion.
KARR-seq utilizes chemical crosslinkers to capture physically proximal transcripts independent
of local RBP concentration and RBP-RNA affinities. 

# Software Pre-requisites

## Dependencies and Prerequisites
The following software are required to run the Snakemake pipeline
+ Python 3.7
+ Cigar (PyPI)
+ Snakemake
+ STAR (>=2.5.2a)
+ Samtools (>=1.1)
+ SeqPrep
+ GNU sort (with --parallel option)

# General instructions to setup pipeline

## 1. Point to the correct executables
Edit the Snakefile by pointing to the correct path of the following variables:
- `STAR_INDEX`: Directory of the STAR index of your target reference FASTA file 
- `S_STAR`
- `S_SAMTOOLS`
- `S_PYTHON`
- `S_SEQPREP`
- `S_SORT`
- `S_PAIRIX`
- `S_BGZIP`

If you are using modules in a cluster environment, you could prepend the module statement prior to the executable location, for e.g. `S_SAMTOOLS = "module samtools/1.1; samtools"`

## 2. Naming convention of FASTQ
Place FASTQ files under the `data/fastq` folder and name them as
`<example1>_P1.fastq.gz` and `<example1>_P2.fastq.gz`

## 3. Run the pipeline
Using the snakemake executable, simply run:
`snakemake -j <no_of_jobs> -s <Snakefile> --latency-wait 120 targets_all`

## 4. Output Files
- Chimeric reads (both gapped and chiastic reads) tab-delimited text file
    - `<dir>/data/pairs/<genome>/MAPQ1_SPAN0/<example1>.dedup.txt.gz`
    - readID chr1 start1 chr2 start2 strand1 strand2 mapq1 mapq2 armlen1 armlen2 order1 order2 ignore
- For fast, random retrieval, use the [Pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) file format from the 4D Nucleome Omics Data Standards Working Group
    - `<dir>/data/pairs/<genome>/MAPQ1_SPAN0/<example1>.dedup.pairs.gz`
- Pairs index (reindex with pairix if outdated)
    - `<dir>/data/pairs/<genome>/MAPQ1_SPAN0/<example1>.dedup.pairs.gz.px2`


# Usage
See the jupyter notebook or HTML rendered notebook for examples on how to utilize the pairs files to cluster chimeric reads, plot contactmaps, arcbands as well as differential contact maps.