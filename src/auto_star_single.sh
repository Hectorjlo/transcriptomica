#!/usr/bin/env bash

# Runs STAR single_end in parallel for each read, hardcoded code.

mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*_1.fastq" | sort)

# Count how many read_files are
n_files=$(( ${#files[@]} ))

# Hardcoded function
run_star() {
    # Arguments
    #   $1: input FASTQ file
    #   $2: output directory prefix (STAR uses a prefix, not a single SAM file)
    #   $3: summary/log directory (unused directly; STAR auto-generates logs)
    STAR --runMode alignReads \
         --genomeDir /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.star/ \
         --readFilesIn "$1" \
         --outSAMtype SAM \
         --outFileNamePrefix "$2" \
         --outSAMunmapped None \
         --runThreadN 2 \
         --outFilterMultimapNmax 7
}
export -f run_star

out_dir="${2%/}/"
for (( i=0; i<${#files[@]}; i+=1 )); do
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%.fastq}"
echo "${files[$i]} ${out_dir}output_${core_name}_ ${out_dir}summary_of_${core_name}.txt"
done | parallel -j ${n_files} --colsep ' ' run_star