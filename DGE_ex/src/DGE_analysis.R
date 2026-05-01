library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)


gene_counts = read.delim("DGE_ex/data/gene-counts.tsv", row.names=1)
annotation = read.delim("DGE_ex/data/annotation.tsv", row.names=1)
gene_name_map = read.delim("DGE_ex/data/mm39-gencode-M36-gene_id-gene_name.txt", header=FALSE, row.names=1)


# Generate factors
# Levels hace que la comparación sea como base m4, entonces todas las
# condiciones contra m4 (male age 4)
age = factor(c(rep("m4",10),rep("m18",10)), levels=c("m4","m18")) 
sex = factor(c(rep("male",5), rep("female",5),rep("male",5), rep("female",5)))
sample_color = c(rep("lightblue",5), rep("lightpink",5),rep("blue",5), rep("red",5))
sample_names = colnames(gene_counts)

# Generate metadata table
meta_data = data.frame(sample_names, age, sex)
meta_data = meta_data %>% remove_rownames %>% column_to_rownames(var="sample_names") # Delete column sample_names and make it the rowname column

# Check that all names in meta_data are in the same order than the columns of gene_counts
all(colnames(gene_counts) == rownames(meta_data))

dds = DESeqDataSetFromMatrix(countData = round(gene_counts), colData = meta_data, design= ~ 0 + age)
# Al hacer el diseño ~ 0 + age, no se automaticamente como base ninguna
# condición, haciendo efectivo lo que se hizo con levels antes
design = model.matrix(~ 0 + age)

dds = DESeqDataSetFromMatrix(countData = round(gene_counts), colData = meta_data, design= ~ 0 + age)
design = model.matrix(~ 0 + age)

vsd = vst(dds)
plotPCA(vsd, intgroup = c("age")) + aes(shape=sex) + theme_classic(base_size=18, base_line_size = 1)

dds = DESeq(dds)

# Add gene length
mcols(dds)$basepairs = annotation[rownames(dds),]

log2_fpkm = log2(fpkm(dds)+0.1) # Calculate log2 FPKM
fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) } # function to convert FPKM to TPM
log2_tpm = apply(log2_fpkm, 2, fpkm2tpm_log2) # We apply the function

# Save nice table
gene_names = gene_name_map[rownames(log2_tpm),]
write.table(cbind(gene_names, log2_tpm), "DGE_ex/results/TPM_log2-table.txt", sep="\t", quote=FALSE)


resultsNames(dds)

contrasts = makeContrasts(m18_vs_m4=agem18-agem4, 
                         levels=design)


res = results(dds, contrast=contrasts[,"m18_vs_m4"])
res$Gene_name = gene_name_map[rownames(res),] # Add gene name

# Select LFC and FDR tresholds
FDR = 1e-10
LFC = 0.5

# Get up, and down regulated genes:

up = (res$log2FoldChange > LFC) & (res$padj < FDR) 
up[which(is.na(up))] = FALSE # eliminate NAs
cat ("Upregulated: ", sum(up), "\n")

down = (res$log2FoldChange < -LFC) & (res$padj < FDR) 
down[which(is.na(down))] = FALSE # eliminate NAs
cat ("Downregulated: ", sum(down), "\n")

# Save up and down tables
write.table(res[up,], paste("DGE_ex/results/deseq-DEG-Up-", FDR, ".txt", sep=""), sep="\t", quote=FALSE, row.names=T)
write.table(res[down,], paste("DGE_ex/results/deseq-DEG-Down-", FDR, ".txt", sep=""), sep="\t", quote=FALSE, row.names=T)


# Asign LFC and FDR limits
volcano_LFC_limit=12
volcano_FDR_limit=70

# Asignar colors
vpcolors = c("gray", "blue", "red") 
names(vpcolors) = c("NO", "DOWN", "UP")

# Add DE column
res$DE="NO" 
res[up,"DE"]="UP"
res[down,"DE"]="DOWN"

p = ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), col=DE)) +
        geom_point() +
        labs(title = "Volcano plot", x = "log2 Expression fold change", y = "-log10 FDR") +
        xlim(-volcano_LFC_limit, volcano_LFC_limit) + ylim(0, volcano_FDR_limit) +
        scale_color_manual(values=vpcolors) +
        geom_vline(xintercept=c(-LFC,LFC), col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(FDR), col="black", linetype = "longdash") +
        theme_classic(base_size = 15, base_line_size = 1)
p



significant = res[up | down,] # Get significant genes
significant_order = significant[order(significant$padj),] # Order by FDR
significant = head(rownames(significant_order), n = 2000) # Get top 2,000
zscore_significant = t(scale(t(log2_tpm[significant,])))  # Calculate z score from TPM
p = Heatmap(zscore_significant, cluster_rows = T, cluster_columns = F, show_row_names = F, name = "Z score", km = 2, column_title = "Heatmap all")
draw(p)


significant = head(rownames(significant_order), n = 20)
zscore_top = t(scale(t(log2_tpm[significant,]))) 
p = Heatmap(zscore_top, cluster_rows = T, cluster_columns = F, row_labels = gene_name_map[rownames(zscore_top),], name = "Z score", km = 2, column_title = "Top 20 significant genes")
draw(p)




assay(vsd) = limma::removeBatchEffect(assay(vsd), batch=sex, design=model.matrix(~ 0 + age))
plotPCA(vsd, intgroup = "age") + aes(shape=sex) + theme_classic(base_size=18, base_line_size = 1)

design(dds) <- ~ 0 + age + sex
design = model.matrix(~ 0 + age + sex)
dds = DESeq(dds)

resultsNames(dds)
contrasts = makeContrasts(m18_vs_m4=agem18-agem4, 
                         levels=design)



res = results(dds, contrast=contrasts[,"m18_vs_m4"])

# Select LFC and FDR tresholds
FDR = 1e-10
LFC = 0.5

# Get up, and down regulated genes:

up = (res$log2FoldChange > LFC) & (res$padj < FDR) 
up[which(is.na(up))] = FALSE # eliminate NAs
cat ("Upregulated: ", sum(up), "\n")

down = (res$log2FoldChange < -LFC) & (res$padj < FDR) 
down[which(is.na(down))] = FALSE # eliminate NAs
cat ("Downregulated: ", sum(down), "\n")
