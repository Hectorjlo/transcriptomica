#!/bin/bash -ue
mkdir -p SRR9127308_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9127308_1.fastq SRR9127308_2.fastq -o SRR9127308_raw_fqc -t 20
