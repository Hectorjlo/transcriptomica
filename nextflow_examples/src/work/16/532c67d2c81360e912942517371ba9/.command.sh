#!/bin/bash -ue
mkdir -p SRR9126859_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9126859_1.fastq clean_SRR9126859_2.fastq -o SRR9126859_cleaned_fqc -t 20
