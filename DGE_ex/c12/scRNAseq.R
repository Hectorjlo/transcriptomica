library(Seurat)
library(ggplot2)
library(celldex)
library(SingleR)
library(dplyr)
library(presto)
# Cargar archivos
matrix_samp8 = ReadMtx(
    mtx = "DGE_ex/c12/Clase12-comandos/single_cell/matrix/matrix.mtx.gz",
    features = "DGE_ex/c12/Clase12-comandos/single_cell/matrix/features.tsv.gz",
    cells = "DGE_ex/c12/Clase12-comandos/single_cell/matrix/barcodes.tsv.gz"
)

# Objeto de Seurat
min_cells = 3 # Include genes with expression in at least 3 cells
min_features = 200 # Include cells with at least 200 detected genes

samp8 = CreateSeuratObject(
    counts = matrix_samp8,
    project = "samp8",
    min.cells = min_cells, 
    min.features = min_features 
)

# Filtrar células muertas y potencionales dobletes
samp8[["mt_percent"]] = PercentageFeatureSet(samp8, pattern = "^mt-")

# To check metadata:
head(samp8@meta.data)


# Ahora, necesitamos decidir dónde ponemos el filtro de los genes con 
# alto contenido mitocondrial. A su vez, necesitamos eliminar las células 
# con demasiadas counts y demasiadas features porque posiblemente son dobletes.
samp8[["all"]] = 1
VlnPlot(samp8, features = c("mt_percent", "nCount_RNA", "nFeature_RNA"), layer = "counts", group.by = "all",  ncol = 3, pt.size = 0)

# Los scatters plots suelen ser más informativos para decidir visualmente un umbral
mt_percent_lim = 5
count_lim = 12000
feature_lim = 4000

plot1 = FeatureScatter(samp8, feature1 = "nCount_RNA", feature2 = "mt_percent", group.by = "all") +
        geom_hline(aes(yintercept = mt_percent_lim, col = "blue"))
plot2 = FeatureScatter(samp8, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "all") +
        geom_hline(aes(yintercept = feature_lim, col = "blue")) + geom_vline(aes(xintercept = count_lim, col="blue"))

plot1 + plot2 


# Aplicamos filtros
samp8
samp8 = subset(samp8, subset = nCount_RNA < count_lim & nFeature_RNA < feature_lim & mt_percent < mt_percent_lim)
samp8

# Encontrar features variables
# Lo siguiente es buscar cuáles son los features más variables. 
# Estas variable features son usadas para encontrar los clusters de poblaciones.

samp8 = NormalizeData(samp8, normalization.method = "LogNormalize", scale.factor = 10000) # Normalize using default values

samp8 = FindVariableFeatures(samp8, selection.method = "vst", nfeatures = 2000) # Find variable genes. These are useful because usually these are genes that change between cell types.

# Ploteo
top20 = head(VariableFeatures(samp8), 20)
plot1 = VariableFeaturePlot(samp8)
plot1 = LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1

# Escalar datos
# Con esta función generamos una métrica escalada por z-score que
#  nos puede ser útil para comparar entre datasets. El PCA se genera 
# a partir de los datos escalados.
all.genes = rownames(samp8)
samp8 = ScaleData(samp8, features = all.genes) # Scaled data is saved in samp8[["RNA"]]$scale.data


# Realizamos un PCA
samp8 = RunPCA(samp8, features = VariableFeatures(object = samp8))
DimPlot(samp8, reduction = "pca") + NoLegend() # full PCA

# Si estamos trabajando con dos sets de datos distintos, 
# este es un buen momento para determinar si necesitas armonizarlos 
# (chequen el paquete harmony). Si ven un efecto obvio de batch effect, 
# este es el momento de integrarlos. Si no lo ven, procedan y 
# vuelvan a checar cómo se ven las cosas en el UMAP.

# Adicionalmente, podemos plotear a los genes que más pesan en cada uno 
# de los loadings del PCA:

VizDimLoadings(samp8, dims = 1:2, reduction = "pca") # top genes for PC1 and PC2

# Los componentes del PCA se utilizan para hacer los clústers de los genes. 
# En este momento deben decidir cuántos componentes usan. Para eso, siempre
#  es útil hacer un elbow plot:
ElbowPlot(samp8) # Make elbow plot (similar to a scree plot, but with points) to see ideal number of clusters

# Usando el elbow plot pueden elegir el número de PCs que contenga la 
# mayor parte de la varianza. Nosotros usaremos 15. Cuando tengan duda, 
# siempre es mejor empezar de más a menos componentes principales.

# Generar clusters y UMAP
res = 0.5 # This is the resolution for finding clusters. 0.5 is adviced by seurat as a starting point

samp8 = FindNeighbors(samp8, dims = 1:15) # Find neighbors based on 15 first PCs
samp8 = FindClusters(samp8, resolution = res)


# El número de clústers que encuentran está determinado principalmente por 
# la resolución. Seurat recomienda comenzar por 0.5, pero no hay una regla específica. 
# Es por eso que el número de clústers es, esencialmente, un dato muy poco informativo.

# Y ahora sí, calculamos el UMAP. En este paso, si lo preferimos, podemos 
# también calcular el TSNE usando la función RunTSNE
samp8 = RunUMAP(samp8, dims = 1:15)
DimPlot(samp8, reduction = "umap")

# Encontrar marcadores
# El concepto de “Análisis de expresión diferencial” es distinto en un 
# bulk y un single-cell. En un single-cell por lo general no decimos que 
# estamos haciendo expresión diferencial, sino que estamos buscando marcadores. 
# Un análisis obvio que debemos realizar es buscar marcadores para todas 
# nuestras poblaciones para que podamos anotarlas. En este caso, nos quedamos sólo
# con los marcadores positivos, es decir, los que están sobre-expresados
# (hay que recordar que por su naturaleza, el scRNA-seq es mucho más
# sensible encontrando sobre-expresión que sub-expresión).

# find markers for every cluster compared to all remaining cells, report only the positive ones

samp8_markers = FindAllMarkers(samp8, only.pos = TRUE)

samp8_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# Con la lista de todos nuestros marcadores podemos hacer un 
# heatmap para ver cómo se expresan en nuestros clústers.

samp8_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 4) %>%
    ungroup() -> top

DoHeatmap(samp8, features = top$gene) + NoLegend()

# FindMarkers es la función fundamental para hacer expresión diferencial. 
# Podríamos por ejemplo, comparar el clúster 5 contra los clústers 0 y 3. 
# Aunque eso lo vamos a retomar más adelante cuando tengamos anotadas 
# nuestras células, este es el comando básico:
cluster5_markers = FindMarkers(samp8, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_markers)


# Anotación
# Hay dos formas de hacer la anotación. Una es manual, 
# analizando los marcadores de cada clúster y determinando qué tipo celular es. 
# Otra (para mi gusto mejor en la mayoría de los casos) es hacerlo de forma
#  automática utilizando un dataset ya anotado.

# Veamos cómo haríamos una anotación manual:
manual_annotation = c("Perro", "Gato", "Canario", "Gorila", "Bonobo", "Tejon", "Humano", "Serpiente", "Leon", "Aguila", "Gorila", "Gato")

names(manual_annotation) = levels(samp8)
samp8 = RenameIdents(samp8, manual_annotation)
DimPlot(samp8, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# As we don´t want this mock annotation, we go back to seurat clusters:
samp8 = FindClusters(samp8, resolution = res)

# Ahora, haremos una anotación automática. Existen diferentes datasets y 
# formas de hacerlo. Una muy sencilla es usar SingleR que tienen una referencia
# basada en distintos datasets de RNAseq de ratón:
ref = fetchReference("mouse_rnaseq", "2024-02-26")
sce = as.SingleCellExperiment(DietSeurat(samp8)) # Seurat object needs to be changed to a sce object

annotation_mouse_rnaseq = SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.main)
rm (sce, ref) # delete sce and ref objects
samp8[["annotation_mouse_rnaseq"]] = annotation_mouse_rnaseq$pruned.labels

DimPlot(samp8, reduction = "umap", group.by = "annotation_mouse_rnaseq")

# Si vemos tipos celulares extraños, podemos crear un filtro para las anotaciones 
# con muy pocas células (menos del 1% de la población total)

t = table(samp8@meta.data$annotation_mouse_rnaseq)
valid_cell_types = rownames(t)[t > 50] # We get the list of samples and cell types with more than 50 cells to subset
keep = samp8$annotation_mouse_rnaseq %in% valid_cell_types
keep_cells = rownames(samp8@meta.data)[keep]

DimPlot(samp8, reduction = "umap", group.by = "annotation_mouse_rnaseq", cells = keep_cells)


# Marcadores por tipo celular
# Ahora que tenemos anotadas las células, podemos recalcular los marcadores, 
# pero esta vez por cada tipo celular
samp8_markers = FindAllMarkers(samp8, only.pos = TRUE, group.by = "annotation_mouse_rnaseq", min.cells.group = 50)
samp8_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

samp8_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 4) %>%
    ungroup() -> top

DoHeatmap(samp8, features = top$gene) + NoLegend()
DoHeatmap(samp8, features = top$gene, cells = keep_cells, group.by = "annotation_mouse_rnaseq") + NoLegend()

# Plots útiles
# Una vez que tengamos identificados genes particulares que nos interesen, 
# podemos hacer plots interesantes con ellos. Por ejemplo, digamos que queremos 
# ver cómo se expresa los marcadores de microglía Hexb y Cx3cr1. Podemos hacer un 
# violinplot:

VlnPlot(subset(samp8, cells = keep_cells), features = c("Hexb", "Cx3cr1"), group.by = "annotation_mouse_rnaseq")

# Ahora supongamos que queremos plotear la expresión en el UMAP:
FeaturePlot(samp8, features = c("Hexb", "Cx3cr1"), cells = keep_cells)

# O si queremos hacer un heatmap:
DoHeatmap(samp8, features = c("Hexb", "Cx3cr1"), cells = keep_cells, group.by = "annotation_mouse_rnaseq")


# DEG
# Es muy posible que ahora que tengamos re-anotadas las células queramos 
# hacer un análisis de expresión diferencial mucho más preciso. Para eso usamos 
# la función FindMarkers. Sólo que ahora probaremos un par de métodos al comparar 
# microglía contra neuronas. Pero primero, vamos a crear una función para hacer 
# volcano plots:
# Function for volcano plot:
DE_volcano <- function(DEG, FDR_lim, LFC_lim, FDR_plot_lim, LFC_plot_lim, table_folder, plot_folder, basename) {
    # Set colors
    vpcolors = c("gray", "blue", "red") 
    names(vpcolors) = c("NO", "DOWN", "UP")
    # Get up and downregulated genes
    up = (DEG$avg_log2FC > LFC_lim) & (DEG$p_val_adj < FDR_lim) # upregulated genes   
    cat ("Upregulated in", basename, sum(up), "\n")
    down = (DEG$avg_log2FC < -LFC_lim) & (DEG$p_val_adj < FDR_lim) # downregulated genes
    cat ("Downregulated in", basename, sum(down), "\n")
    
    # Create DE column for filtering and color
    DEG$DE="NO" 
    DEG[up,"DE"]="UP"
    DEG[down,"DE"]="DOWN"

    # Write tables
    write.table(DEG[up,], paste(table_folder, basename, "-Up.txt", sep=""), sep="\t", quote=FALSE, row.names=T, col.names=T)
    write.table(DEG[down,], paste(table_folder, basename, "-Down.txt", sep=""), sep="\t", quote=FALSE, row.names=T, col.names=T)

    # Make volcano plot
    ggplot(data=DEG, aes(x=avg_log2FC, y=-log10(p_val_adj), col=DE)) +
        geom_point() +
        labs(title = paste("Volcano plot ", basename, sep=""), x = "log2 Expression fold change", y = "-log10 FDR") +
        xlim(-LFC_plot_lim, LFC_plot_lim) + ylim(0, FDR_plot_lim) +
        scale_color_manual(values=vpcolors) +
        geom_vline(xintercept=c(-LFC_lim,LFC_lim), col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(FDR_lim), col="black", linetype = "longdash") +
        theme_classic(base_size = 15, base_line_size = 1)
    ggsave(paste(plot_folder, basename, "-volcano_plot.pdf", sep=""))
}

# Ahora probamos la expresión diferencial usando un simple test de wilcoxon que es el default:

FDR = 1e-5
LFC = 0.5
volcano_LFC_limit=10
volcano_FDR_limit=300

DEG = FindMarkers(samp8, ident.1 = "Microglia", ident.2 = "Neurons", group.by= "annotation_mouse_rnaseq", test.use = "wilcox") # Note that ident.1 is the numerator and ident.2 is the denominator in the comparison. So, when something is upregulated, it is upregulated in ident.1 compared to ident.2
DE_volcano(DEG, FDR, LFC, volcano_FDR_limit, volcano_LFC_limit, "./", "./", "wilcox")

# Sin embargo, es mucho más conveniente usar algoritmos más avanzados. 
# MAST y DESeq2 están disponibles para usarse de manera fácil y rápida. 
# De los dos, MAST es el más rápido y sensible.

# Wilcoxon
DEG = FindMarkers(samp8, ident.1 = "Microglia", ident.2 = "Neurons", group.by= "annotation_mouse_rnaseq", test.use = "MAST")
DE_volcano(DEG, FDR, LFC, volcano_FDR_limit, volcano_LFC_limit, "./", "./", "MAST")

# DESeq2
DEG = FindMarkers(samp8, ident.1 = "Microglia", ident.2 = "Neurons", group.by= "annotation_mouse_rnaseq", test.use = "DESeq2")
DE_volcano(DEG, FDR, LFC, volcano_FDR_limit, volcano_LFC_limit, "./", "./", "DESeq2")

# Últimos consejos
# Hay mucho más que aprender del análisis de single cell, 
# desafortunadamente hay poco tiempo. El campo evoluciona muy rápido y la 
# sintaxis de seurat a veces es inconsistente. No les importa romper cosas 
# con tal de mejorar. Así que si algo ya no funciona, probablemente es un 
# parámetro que ya se movió de nombre o de lugar. Afortunadamente, encontrar 
# la solución es más fácil que nunca.
