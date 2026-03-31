#!/usr/bin/env bash

# Runs hisat2 paired_end in parallel for each read, harcoded code.

mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*.fastq" | sort)

# Count how many read_files are / 2, pe reads
n_files=$(( ${#files[@]} / 2 ))


# Hardcoded function
run_hisat2() {
    # Arguments
    #   -x: 
    #   -U:
    #   -S:
    #   --summary-file:
    #   -k:
    #   --no-unal:
    #   --no-mixed:
    hisat2 -x /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat \
           -1 "$1" \
           -2 "$2" \
           -S "$3" \
           --summary-file "$4" \
           -k 7 \
           --no-unal \
           --no-mixed \
           --no-discordant

}
export -f run_hisat2

# For each read_file do ... until i is less than the number of elements in the array "files"
# starting with i = 0
out_dir="${2%/}/"
for (( i=0; i<${#files[@]}; i+=2 )); do
# The echo will send 3 strings:
#   1st -> The fastp file of unpaired reads 
#   2nd -> Name of the output SAM file
#   3rd -> Name of the summary file 
# These strings are build in the echo
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%_1.fastq}"
echo "${files[$i]} ${files[$i+1]} ${out_dir}output_${core_name}_pe.sam ${out_dir}summary_of_${core_name}_pe.txt"
done | parallel -j ${n_files} --colsep ' ' run_hisat2