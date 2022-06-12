# AppliedPlantGenomics (GE31/MM12)
This repository contains materials and scripts used for Applied Plant Genomics.

## Calculate Sequence Read Stats ##
This script calculates several statistics of a given FASTQ file: number of reads, number of sequenced nucleotides, average read length (mean), GC content, and N50 length of the reads. The script can handle uncompressed and compressed FASTQ files.


python3 FASTQ_stats3.py
--in <FASTQ_FILE>

## _De novo_ Genome Sequence Assembly ##
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


## Run AUGUSTUS for Gene Prediction ###
[AUGUSTUS]() is one of the standard tools for the generation of _ab initio_ gene predictions. 


augustus \
--gff3=on  \
--UTR=on \
--uniqueGeneId=true \
--codingseq=on\
--species=arabidopsis \
ASSEMBLY \
> annotation.gff

## Extract Peptide and Coding Sequences ##
There is an [AUGUSTUS-associated perl script](https://bioinf.uni-greifswald.de/augustus/binaries/scripts/) for the extraction of peptide and coding sequences called getAnnoFasta.pl. This script requires the genome sequence assembly file (FASTA) and the annotation file (GFF) as input.

getAnnoFasta.pl \
--seqfile=ASSEMBLY_FILE \
ANNOTATION_GFF_FILE


## Add Functional Annotations to Predicted Genes ##
This script allows to assign functional annotation terms to predicted genes based on sequence similarity to previously characterized sequences.

python3 construct_anno3.py \
--out OUTPUT_FOLDER \
--in PREDICTED_PEPTIDE_FILE \
--ref ath.pep.fasta \
--anno ath.anno.txt


## Long Read Mapping ##
[minimap2]() can be used to align long reads to a genome sequence.


/vol/data/tools/minimap2-2.24_x64-linux/minimap2 \
-ax map-ont --secondary=no -t 10 \
GENOMSE_ASSEMBLY_FASTA \
FASTQ_FILE \
> MAPPING_FILE


samtools view -Sb \
SAM_FILE \
> BAM_FILE


samtools sort -@ 10 \
-o OUTPUT_BAM_FILE \
INPUT_BAM_FILE


samtools index <BAM_FILE>



## References ##



