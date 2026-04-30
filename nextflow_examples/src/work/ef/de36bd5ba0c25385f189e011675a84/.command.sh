#!/bin/bash -ue
mkdir -p SRR9126934_cleaned_fqc
/export/apps/bioconda/bin/fastqc clean_SRR9126934_1.fastq clean_SRR9126934_2.fastq -o SRR9126934_cleaned_fqc -t 20
