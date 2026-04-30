#!/bin/bash -ue
mkdir -p SRR9127293_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9127293_1.fastq SRR9127293_2.fastq -o SRR9127293_raw_fqc -t 20
