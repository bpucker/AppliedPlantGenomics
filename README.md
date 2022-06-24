# Applied Plant Genomics (GE31/MM12)
This repository contains materials and scripts used for Applied Plant Genomics.

## Calculate Sequence Read Stats ##
This script calculates several statistics of a given FASTQ file: number of reads, number of sequenced nucleotides, average read length (mean), GC content, and N50 length of the reads. The script can handle uncompressed and compressed FASTQ files.



```
Usage:
  python3 FASTQ_stats3.py
 
Mandatory:
  --in   STR   FASTQ file
```



## _De novo_ Genome Sequence Assembly ##
[Canu]() is an established tool for the generation of long read plant genome sequence assemblies. Alternatives are miniasm, Flye, and shasta.

Here is an explanation how to run Canu:


```
/vol/data/tools/canu-2.2/bin/canu \
-d /vol/data/members/<USER_NAME>/assembly/ \
-p <ASSEMBLY_NAME> \
genomeSize=1m \
-nanopore /vol/data/data/reads.fastq.gz \
minReadLength=10000 \
> /vol/data/members/<USER_NAME>/docu.txt \
2> /vol/data/members/<USER_NAME>/err.txt &
```

## Calculate Assembly Statistics ##
This script calculates general assembly statistics and removes short contigs below a certain cutoff.

```
Usage:
  python3 contig_stats3.py
 
Mandatory:
  --input           STR    Assemblyfile
 		
Optional:
  --out              STR    Output folder
  --min_contig_len   INT   Minimal contig length [500]
```


## Run AUGUSTUS for Gene Prediction ###
[AUGUSTUS]() is one of the standard tools for the generation of _ab initio_ gene predictions. 

```
augustus \
--gff3=on  \
--UTR=on \
--uniqueGeneId=true \
--codingseq=on \
--species=arabidopsis \
ASSEMBLY \
> annotation.gff
```


## Extract Peptide and Coding Sequences ##
There is an [AUGUSTUS-associated perl script](https://bioinf.uni-greifswald.de/augustus/binaries/scripts/) for the extraction of peptide and coding sequences called getAnnoFasta.pl. This script requires the genome sequence assembly file (FASTA) and the annotation file (GFF) as input.

```
getAnnoFasta.pl \
--seqfile=ASSEMBLY_FILE \
ANNOTATION_GFF_FILE
```

## Add Functional Annotations to Predicted Genes ##
This script allows to assign functional annotation terms to predicted genes based on sequence similarity to previously characterized sequences.

```
Usage:
  python3 construct_anno3.py
 
Mandatory:
  --out OUTPUT_FOLDER \
  --in PREDICTED_PEPTIDE_FILE \
  --ref ath.pep.fasta \
  --anno ath.anno.txt
```


## Long Read Mapping ##
[minimap2](https://github.com/lh3/minimap2) can be used to align long reads to a genome sequence.

```
/vol/data/tools/minimap2-2.24_x64-linux/minimap2 \
-ax map-ont --secondary=no -t 10 \
GENOMSE_ASSEMBLY_FASTA \
FASTQ_FILE \
> MAPPING_FILE
```

```
samtools view -Sb \
SAM_FILE \
> BAM_FILE
```

```
samtools sort -@ 10 \
-o OUTPUT_BAM_FILE \
INPUT_BAM_FILE
```

```
samtools index <BAM_FILE>
```


## References ##

Pucker B, Holtgräwe D, Rosleff Sörensen T, Stracke R, Viehöver P, Weisshaar B (2016) A De Novo Genome Sequence Assembly of the Arabidopsis thaliana Accession Niederzenz-1 Displays Presence/Absence Variation and Strong Synteny. PLoS ONE 11(10): e0164321. doi:[10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321).

Pucker B, Holtgräwe D, Stadermann KB, Frey K, Huettel B, Reinhardt R, et al. (2019) A chromosome-level sequence assembly reveals the structure of the Arabidopsis thaliana Nd-1 genome and its gene set. PLoS ONE 14(5): e0216233. doi:[10.1371/journal.pone.0216233](https://doi.org/10.1371/journal.pone.0216233).

Pucker, B., Kleinbölting, N. & Weisshaar, B. Large scale genomic rearrangements in selected Arabidopsis thaliana T-DNA lines are caused by T-DNA insertion mutagenesis. BMC Genomics 22, 599 (2021). doi:[10.1186/s12864-021-07877-8](https://doi.org/10.1186/s12864-021-07877-8).

