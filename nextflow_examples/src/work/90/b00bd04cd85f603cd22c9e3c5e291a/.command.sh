#!/bin/bash -ue
mkdir -p SRR9127396_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9127396_1.fastq clean_SRR9127396_2.fastq -o SRR9127396_cleaned_fqc -t 20
