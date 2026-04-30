#!/bin/bash -ue
mkdir -p SRR9126859_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9126859_1.fastq SRR9126859_2.fastq -o SRR9126859_raw_fqc -t 20
