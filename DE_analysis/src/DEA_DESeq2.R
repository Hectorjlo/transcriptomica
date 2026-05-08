# Se cargaron las librerías necesarias para el análisis
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)

# Se obtuvieron los argumentos del paseo del CLI
# args <- commandArgs(trailingOnly = TRUE)
# gene_counts_path <- args[[1]]
# annotacion_file_path <- args[[2]]
# gene_name_map_file_path <- args[[3]]

# Fijo ejemplo para hisat/single_end
gene_counts_path <- "DE_analysis/results/hisat2/feature_counts/single_end/feature_counts_global.tsv"
annotacion_file_path <- ""
gene_name_map_file_path <- ""


# Leer los archivos de conteos, anotación y del genemap
gene_counts = read.delim(gene_counts_path, row.names=1)
annotation = read.delim(annotacion_file_path, row.names=1)
gene_name_map = read.delim(gene_name_map_file_path, header=FALSE, row.names=1)