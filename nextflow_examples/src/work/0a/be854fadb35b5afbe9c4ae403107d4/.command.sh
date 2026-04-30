#!/bin/bash -ue
mkdir -p SRR9127396_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9127396_1.fastq SRR9127396_2.fastq -o SRR9127396_raw_fqc -t 20
