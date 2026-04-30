#!/bin/bash -ue
/export/apps/bioconda/envs/fastp/bin/fastp -i SRR9127396_1.fastq -I SRR9127396_2.fastq           -o clean_SRR9127396_1.fastq -O clean_SRR9127396_2.fastq           --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20           --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20           --detect_adapter_for_pe --trim_poly_g -l 50 -w 20
