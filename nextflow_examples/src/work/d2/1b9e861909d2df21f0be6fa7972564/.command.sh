#!/bin/bash -ue
mkdir -p SRR9126934_raw_fqc
/export/apps/bioconda/bin/fastqc SRR9126934_1.fastq SRR9126934_2.fastq -o SRR9126934_raw_fqc -t 20
