#!/bin/bash -ue
mkdir -p SRR9127063_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9127063_1.fastq clean_SRR9127063_2.fastq -o SRR9127063_cleaned_fqc -t 20
