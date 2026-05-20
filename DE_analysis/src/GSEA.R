library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)


res <- readRDS(file = "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/results.rds")

res <- res[order(-res$stat), ]

gene_list <- res$stat
names(gene_list) <- rownames(res)

gse <- gseGO(geneList = gene_list,
            ont = "BP",
            keyType = "ENSEMBL",
            OrgDb = "org.Mm.eg.db",
            eps = 1e-300
        )

saveRDS(gse, file = "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/gse.rds")

p <- gseaplot2(gse,
               geneSetID = c(grep("osteoclast proliferation", as.data.frame(gse)$Description), 
                            grep("negative regulation of osteoclast differentiation", as.data.frame(gse)$Description)),
               title = "",
               base_size = 14,
               pvalue_table = TRUE,
               pvalue_table_columns = c("NES", "p.adjust")
            )

p[[1]] <- p[[1]] +  geom_hline(yintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.5)

p

ggsave(filename = "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/plots/GSEA_plot.png",
    plot = p, 
    dpi = 1200, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)
