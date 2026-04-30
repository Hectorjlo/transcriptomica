#!/bin/bash -ue
mkdir -p SRR9127357_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9127357_1.fastq SRR9127357_2.fastq -o SRR9127357_raw_fqc -t 20
