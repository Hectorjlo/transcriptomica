#!/bin/bash -ue
mkdir -p SRR9127308_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9127308_1.fastq clean_SRR9127308_2.fastq -o SRR9127308_cleaned_fqc -t 20
