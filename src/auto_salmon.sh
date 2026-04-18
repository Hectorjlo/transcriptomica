#!/usr/bin/env bash
# Runs salmon quant for single and paired reads
# This script expects three arguments
# Usage:
#   ./auto_salmon.sh <input_dir> <output_dir_se> <output_dir_pe>

# Create an array of the names that match: *.sam
mapfile -t single < <(find "$1" -maxdepth 1 -type f -name "clean_*_1.fastq" | sort) # Single_end
mapfile -t paired < <(find "$1" -maxdepth 1 -type f -name "clean_*.fastq" | sort) # Paired_end
# Calculate total file length
n_single=$(( ${#single[@]} )) # Single_end
n_paired=$(( ${#paired[@]} / 2 )) # Paired_end

is_single=true
# Hardcoded function
run_salmon() {
    # Arguments
    #   -i: Index file
    #   -l: Mode of Alignment
    #   -r: Read file
    #   -o: Output file path
    if $4; then
        salmon quant -i /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.salmon \
                    -l A \
                    -r "$1" \
                    -o "$2"
    else
        salmon quant -i /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.salmon \
                    -l A \
                    -1 "$1" \
                    -2 "$2" \
                    -o "$3"
    fi


}
# Export the function, needed for parallel to use it
export -f run_salmon

# Take the second argument passed to the main call ($2) and remove a "/" adding at the end a "/"
out_dir="${2%/}/"

for (( i=0; i<${#single[@]}; i+=1 )); do
# Build the core_name var, it takes the element "i" of the array "single"
# starting from the front of the string will eliminate all "*" until a "/" is found,
# then again starting from the front will eliminate all ultil a "clean_" match is found,
# finally starting from the back of the string will eliminate the ".fastq"
core_name="${single[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%.fastq}"
# The echo will send 4 strings:
#   1st -> The input file
#   2nd -> Output name of the file
#   3rd -> Placeholder
#   4th -> Boolean-like flag
# These strings are build in the echo
echo "${single[$i]} ${out_dir}output_${core_name}.salmon _ $is_single"
done | parallel -j ${n_single} --colsep ' ' run_salmon


## Paired
# Set flag to false
is_single=false
# Take the third argument in the main call
out_dir="${3%/}/"
for (( i=0; i<${#paired[@]}; i+=2 )); do

core_name="${paired[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%_1.fastq}"
# The echo will send 4 strings:
#   1st -> The input file read 1 
#   2nd -> The input file read 2
#   3rd -> Output name of the file
#   4th -> Boolean-like flag
# These strings are build in the echo
echo "${paired[$i]} ${paired[$i+1]} ${out_dir}output_${core_name}.salmon $is_single"
done | parallel -j ${n_paired} --colsep ' ' run_salmon

