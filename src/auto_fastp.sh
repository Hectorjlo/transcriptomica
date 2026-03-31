#!/usr/bin/env bash

# Runs fastp in parallel for n-fastq/2 files with hardcoded options

# Using mapfile to create an array of the files
# find and sort are used to keep only files but folders
mapfile -t files < <(find ../SRRs -maxdepth 1 -type f -name "*.fastq" | sort)
#          ^^^^ Name of the array

# To run parallel it should take half of the fastq files
n_pairs=$(( ${#files[@]} / 2 ))
#             ^^^^^^^^^ Accessing the number of elements in the array
            
# Function to execute fastp for a pair
run_fastp() {
    fastp -i "$1" -I "$2" \
          -o clean_$3 -O clean_$4 \
          --adapter_sequence AAGCAGTGGTATCAACGCAGAGTAC \
          --adapter_sequence_r2 AAGCAGTGGTATCAACGCAGAGTAC \
          --trim_front1 15 \
          --trim_front2 15 \
          --detect_adapter_for_pe --trim_poly_g -l 50
}
export -f run_fastp

# For each pair do ... until i is less than the number of elements in the array "files"
# starting with i = 0
for (( i=0; i<${#files[@]}; i+=2 )); do
#                           ^^^ Pass to the next pair
# The echo send 4 strings, the first two are the path of the files fastq
# the next two are only the names of those files (deletes the ../SRRs/ part keeps the name)
    echo "${files[$i]} ${files[$i+1]} ${files[$i]#../SRRs/} ${files[$i+1]#../SRRs/}" 
done | parallel -j "$n_pairs" --colsep ' ' run_fastp  
#      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For each pair (separated by spaces) execute the function
# By this way each fastp is run in a separate theread, concurrency at it most


