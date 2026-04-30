#!/bin/bash -ue
mkdir -p SRR9127293_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9127293_1.fastq clean_SRR9127293_2.fastq -o SRR9127293_cleaned_fqc -t 20
