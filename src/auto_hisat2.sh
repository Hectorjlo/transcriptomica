#!/usr/bin/env bash

# Runs hisat2 single_end in parallel for each read, harcoded code.

mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*_1.fastq" | sort)

# Count how many read_files are
n_files=$(( ${#files[@]} ))


# Hardcoded function
run_hisat2() {
    # Arguments
    #   -x: 
    #   -U:
    #   -S:
    #   --summary-file:
    #   -k:
    #   --no-unal:
    hisat2 -x /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat \
           -U "$1" \
           -S "$2" \
           --summary-file "$3" \
           -k 7 \
           --no-unal

}
export -f run_hisat2

# For each read_file do ... until i is less than the number of elements in the array "files"
# starting with i = 0
out_dir="${2%/}/"
for (( i=0; i<${#files[@]}; i+=1 )); do
# The echo will send 3 strings:
#   1st -> The fastp file of unpaired reads 
#   2nd -> Name of the output SAM file
#   3rd -> Name of the summary file 
# These strings are build in the echo
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%.fastq}"
echo "${files[$i]} ${out_dir}output_${core_name}.sam ${out_dir}summary_of_${core_name}.txt"
done | parallel -j ${n_files} --colsep ' ' run_hisat2