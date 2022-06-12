# AppliedPlantGenomics (GE31/MM12)
This repository contains materials and scripts used for Applied Plant Genomics.

## Calculate Sequence Read Stats ##
This script calculates several statistics of a given FASTQ file: number of reads, number of sequenced nucleotides, average read length (mean), GC content, and N50 length of the reads. The script can handle uncompressed and compressed FASTQ files.


python3 FASTQ_stats3.py
--in <FASTQ_FILE>

## De novo Genome Sequence Assembly ##
[Canu]() is an established tool for the generation of long read plant genome sequence assemblies. Alternatives are miniasm, Flye, and shasta.

Here is an explanation how to run Canu:

/vol/data/tools/canu-2.2/bin/canu \
-d /vol/data/members/<USER_NAME>/assembly/ \
-p <ASSEMBLY_NAME> \
genomeSize=1m \
-nanopore /vol/data/data/reads.fastq.gz \
minReadLength=10000 \
> /vol/data/members/<USER_NAME>/docu.txt \
2> /vol/data/members/<USER_NAME>/err.txt &


## Calculate Assembly Statistics ##
This script calculates general assembly statistics and removes short contigs below a certain cutoff.


python3 contig_stats3.py
--input <ASSEMBLY_FILE>
--out <OUTPUT_FOLDER>
--min_contig_len <MINIMAL_CONTIG_LENGTH>[500]



