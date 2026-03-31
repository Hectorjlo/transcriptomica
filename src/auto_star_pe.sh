#!/usr/bin/env bash

# Runs STAR paired_end in parallel for each read, hardcoded code.

mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*.fastq" | sort)

# Count how many read_files are / 2, pe reads
n_files=$(( ${#files[@]} / 2 ))

# Hardcoded function
run_star() {
    # Arguments
    #   $1: R1 FASTQ file
    #   $2: R2 FASTQ file
    #   $3: output prefix
    #   $4: (unused directly; STAR auto-generates logs)
    STAR --runMode alignReads \
         --genomeDir /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.star/ \
         --readFilesIn "$1" "$2" \
         --outSAMtype SAM \
         --outFileNamePrefix "$3" \
         --outSAMunmapped None \
         --runThreadN 2 \
         --outFilterMultimapNmax 7 \
         --outSAMattributes NH HI AS NM 
}
export -f run_star

out_dir="${2%/}/"
for (( i=0; i<${#files[@]}; i+=2 )); do
    core_name="${files[$i]##*/}"
    core_name="${core_name#clean_}"
    core_name="${core_name%_1.fastq}"
    echo "${files[$i]} ${files[$i+1]} ${out_dir}output_${core_name}_pe_ ${out_dir}summary_of_${core_name}_pe.txt"
done | parallel -j ${n_files} --colsep ' ' run_star