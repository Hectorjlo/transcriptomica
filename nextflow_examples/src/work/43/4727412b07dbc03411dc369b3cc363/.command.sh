#!/bin/bash -ue
mkdir -p SRR9127063_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9127063_1.fastq SRR9127063_2.fastq -o SRR9127063_raw_fqc -t 20
