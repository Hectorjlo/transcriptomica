#!/bin/bash -ue
mkdir -p SRR9127357_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9127357_1.fastq clean_SRR9127357_2.fastq -o SRR9127357_cleaned_fqc -t 20
