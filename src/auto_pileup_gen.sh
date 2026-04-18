#!/usr/bin/env bash

# Runs bamCoverage for paired_end files of expected star and hisat2 alignments
# This script expects four arguments
# Usage:
#   ./auto_pileup_gen.sh <input_star_dir> <input_hisat_dir> <output_star_dir> <output_hisat_dir>
#   ./auto_pileup_gen.sh ../data/trimmed_fastqs ../results/star/paired_end

# Create arrays of the names that match: *.bam files
mapfile -t star_files < <(find "$1" -maxdepth 1 -type f -name "*.bam" | sort) # For Star
mapfile -t hisat_files < <(find "$2" -maxdepth 1 -type f -name "*.bam" | sort) # For Hisat2
# Count how many files where found
n_files_star=$(( ${#star_files[@]} )) # For Star
n_files_hisat=$(( ${#hisat_files[@]} )) # For Hisat2

# Arguments
    #   -p 5: use up to 5 threads per bamCoverage process
    #   -b "$1": input BAM file
    #   -o "$2": output BigWig (.bw) file
    #   -bs 20: bin size of 20 bp
    #   --blackListFileName: ignore genomic regions listed in blacklist BED
    #   --normalizeUsing BPM: normalize coverage to bins per million mapped reads
    #   --skipNAs: skip bins with no coverage data
    #   --ignoreDuplicates: do not count duplicate reads
    #   --samFlagInclude 64: include first mate reads only in paired-end data
run_bamcov_star() {
    bamCoverage -p 5 \
         -b "$1" \
         -o "$2" \
         -bs 20 \
         --blackListFileName /home/hectorjl/4to/transcriptomica/data/exclude_ranges_bed/mm39.excluderanges.bed \
         --normalizeUsing BPM \
         --skipNAs \
         --ignoreDuplicates \
         --samFlagInclude 64
}
run_bamcov_hisat() {
    bamCoverage -p 5 \
         -b "$1" \
         -o "$2" \
         -bs 20 \
         --blackListFileName /home/hectorjl/4to/transcriptomica/data/exclude_ranges_bed/mm39.excluderanges.bed \
         --normalizeUsing BPM \
         --skipNAs \
         --ignoreDuplicates \
         --samFlagInclude 64
}
# Export the function, needed for parallel to use it
export -f run_bamcov_star
export -f run_bamcov_hisat

# Gather file output paths
out_dir_star="${3%/}/" # For Star
out_dir_hisat="${4%/}/" # For Hisat2

# For STAR:
# For each read_file do ... until i is less than the number of elements in the array 
# starting with i = 0, advancing by 1
for (( i=0; i<${#star_files[@]}; i+=1 )); do
# Build the core_name var, it takes the element "i" of the array 
# starting from the front of the string will eliminate all until a "/" is found,
# then again starting from the front will eliminate the first instance of "output_",
# finally starting from the back of the string will eliminate the "_pe_Aligned.out.bam"
core_name="${star_files[$i]##*/}"
core_name="${core_name#output_}"
core_name="${core_name%_pe_Aligned.out.bam}"
# The echo will send 2 strings:
#   1st -> Input bam file
#   2nd -> Path of output file
# These strings are build in the echo
echo "${star_files[$i]} ${out_dir_star}bigwig_${core_name}.bw"
done | parallel -j ${n_files_star} --colsep ' ' run_bamcov_star

# For HISAT2:
for (( i=0; i<${#hisat_files[@]}; i+=1 )); do
# Build the core_name var, it takes the element "i" of the array 
# starting from the front of the string will eliminate all until a "/" is found,
# then again starting from the front will eliminate the first instance of "output_",
# finally starting from the back of the string will eliminate the "_pe.bam"
core_name="${hisat_files[$i]##*/}"
core_name="${core_name#output_}"
core_name="${core_name%_pe.bam}"
echo "${hisat_files[$i]} ${out_dir_hisat}bigwig_${core_name}.bw"
done | parallel -j ${n_files_hisat} --colsep ' ' run_bamcov_hisat