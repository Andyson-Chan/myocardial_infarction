install.packages("Seurat")
require(Seurat)

setwd("/FASTQ_by_lib/MI_nuclei/outs/multi/count/raw_feature_bc_matrix")

#blue print for cellplex multiplexing
#pfade etc. müssen alle adjustiert werden.


data_dir <- '/FASTQ_by_lib/MI_nuclei/outs/multi/count/raw_feature_bc_matrix'

list.files(data_dir)

expression_matrix <- Read10X(data.dir = data_dir)



gd.umis <- expression_matrix$`Gene Expression`

### set UMI cutoff ###
gd.umis <- gd.umis[, colSums(gd.umis)>1000]
######################

gd.htos <- as.matrix(expression_matrix$`Multiplexing Capture`)

gd.htos <- gd.htos[rownames(gd.htos) %in% c("CMO301","CMO302","CMO303", "CMO304", "CMO305"),]



### try removing all cells with no CMO
# gd.htos <- gd.htos [, colSums(gd.htos)>0]
###

joint.bcs <- intersect(colnames(gd.umis), colnames(gd.htos))



gd.umis <- gd.umis[, joint.bcs]

gd.htos <- as.matrix(gd.htos[, joint.bcs])

rownames(gd.htos)



gd.hashtag <- CreateSeuratObject(counts = gd.umis)

gd.hashtag <- NormalizeData(gd.hashtag)

gd.hashtag <- FindVariableFeatures(gd.hashtag)

gd.hashtag <- ScaleData(gd.hashtag, features = VariableFeatures(gd.hashtag))



gd.hashtag[["HTO"]] <- CreateAssayObject(counts = gd.htos)

gd.hashtag <- NormalizeData(gd.hashtag, assay = "HTO", normalization.method = "CLR")


# adjust positive.quantile value to change stringency of doublet/singlet definition
# gd.hashtag <- HTODemux(gd.hashtag, assay = "HTO", positive.quantile = 0.99)
# gd.hashtag <- HTODemux(gd.hashtag, assay = "HTO", positive.quantile = 0.99)
gd.hashtag <- HTODemux(gd.hashtag, assay = "HTO", positive.quantile = 0.99)

table(gd.hashtag$HTO_classification.global)



Idents(gd.hashtag) <- "HTO_maxID"

RidgePlot(gd.hashtag, assay = "HTO", features = rownames(gd.hashtag[["HTO"]])[1:8], ncol = 4)

Idents(gd.hashtag) <- "HTO_classification.global"

VlnPlot(gd.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


# First, we will remove negative cells from the object

gd.hashtag.subset <- subset(gd.hashtag, idents = "Negative", invert = TRUE)

VlnPlot(gd.hashtag.subset, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
hist(gd.hashtag.subset$HTO["CMO305",], breaks = 200)


# Calculate a distance matrix using HTO

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = gd.hashtag.subset, assay = "HTO"))))



# Calculate tSNE embeddings with a distance matrix

gd.hashtag.subset <- RunTSNE(gd.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)

DimPlot(gd.hashtag.subset)
FeaturePlot(gd.hashtag.subset, features = c("Pecam1", "Ptprc", "Pdgfra", "Mcam"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.hashtag.subset, features = c("hto_CMO301", "hto_CMO302", "hto_CMO303", "hto_CMO304", "hto_CMO305"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.hashtag.subset, features = c("Itgax"), min.cutoff = "q9", cols = c("grey90","red4"))
#####################################

##########################################

HTOHeatmap(gd.hashtag, assay = "HTO")

dim(gd.hashtag$HTO)



# Extract the singlets

gd.singlet <- subset(gd.hashtag, idents = "Singlet")

VlnPlot(gd.singlet, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

gd.singlet[["percent.mt"]] <- PercentageFeatureSet(gd.singlet, pattern = "^mt-")

VlnPlot(gd.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = T)



#gd.singlet <- subset(gd.singlet, subset = nCount_RNA < 20000 & percent.mt < 10)

gd.singlet <- NormalizeData(gd.singlet, normalization.method = "LogNormalize", scale.factor = 10000)

gd.singlet <- FindVariableFeatures(gd.singlet, selection.method = "vst", nfeatures = 2000)



# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(gd.singlet), 10)



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(gd.singlet)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1 + plot2



gd.singlet <- ScaleData(gd.singlet)

gd.singlet <- RunPCA(gd.singlet, features = VariableFeatures(object = gd.singlet))

print(gd.singlet[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(gd.singlet, dims = 1:2, reduction = "pca")



gd.singlet <- FindNeighbors(gd.singlet, dims = 1:30)

gd.singlet <- FindClusters(gd.singlet)

gd.singlet <- RunUMAP(gd.singlet, dims = 1:30)

DimPlot(gd.singlet, reduction = "umap", label = T)



set.seed(123)

clust.col = sample(rainbow(25))

DimPlot(gd.singlet, reduction = "umap", label = T,cols = clust.col)

DimPlot(gd.singlet, group.by = "HTO_classification",reduction = "umap",cols = c("red","cyan","green","blue"))



FeaturePlot(gd.singlet, features = c("Pdgfra","Pecam1","Ptprc", "Mcam"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet, features = c("hto_CMO301", "hto_CMO302", "hto_CMO303", "hto_CMO304", "hto_CMO305"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet, features = c("Ttn"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet, features = c("Cd44"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet, features = c("Ctla4", "Mb", "Mdh2", "Atp5b", "Nkx2-5", "Sox4", "Cenpa", "S100a6", "Erbb2", "Tnni1", "Actc1", "Tmsb10", "Zfp706"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet, features = c("Actn2", "Cacna1c", "Kcnj3", "Myh6", "Scn5a", "Tnnt2"), min.cutoff = "q9", cols = c("grey90","red4"))

VlnPlot(gd.singlet, features = c("nCount_RNA"), pt.size = 0.5)
VlnPlot(gd.singlet, features = c("nFeature_RNA"), pt.size = 0.5)


dediff_marker_gene_list <- list(c("Ctla4", "Mb", "Mdh2", "Atp5b", "Nkx2-5", "Sox4", "Cenpa", "S100a6", "Erbb2", "Tnni1", "Actc1", "Tmsb10", "Zfp706"))
dediff_marker_TFgene_list <- list(c("Nkx2-5", "Sox4", "Cenpa", "Tmsb10", "Zfp706"))
cycling_marker_gene_list <- list(c("Mki67", "Ccnd3", "Cdk14", "Knl1", "Kif11", "Cdk14", "Pcna"))
ACM_marker_gene_list <- list(c("Actn2", "Cacna1c", "Kcnj3", "Myh6", "Scn5a", "Tnnt2"))

gd.singlet <- AddModuleScore(object = gd.singlet, features = dediff_marker_gene_list, name = "Dedifferentiation_score")
gd.singlet <- AddModuleScore(object = gd.singlet, features = dediff_marker_TFgene_list, name = "Dedifferentiation_TF_score")
gd.singlet <- AddModuleScore(object = gd.singlet, features = cycling_marker_gene_list, name = "Cycling_score")
gd.singlet <- AddModuleScore(object = gd.singlet, features = ACM_marker_gene_list, name = "ACM_score")

FeaturePlot(object = gd.singlet, features = "ACM_score1")
FeaturePlot(object = gd.singlet, features = "Dedifferentiation_score1")

# Visualize co-expression of two features simultaneously
FeaturePlot(gd.singlet.CM, features = c("Actc1", "Nkx2-5"), blend = TRUE)



# cells in clusters
table(gd.singlet@meta.data$seurat_clusters)

BiocManager::install('limma')
cluster.markers <- FindMarkers(gd.singlet, ident.1 = 7, logfc.threshold = 0.25)

specific.cluster.markers <- FindMarkers(gd.singlet, ident.1 = 0, ident.2 = c(1,6,13,17,18), max.cells.per.ident = 1000)

















#subsetting CM only
length(WhichCells(gd.singlet, idents = c(0, 1, 2, 3, 7)))
gd.singlet.CM <- subset(gd.singlet, idents = c(0,1,2,3,7))

gd.singlet.CM[["percent.mt"]] <- PercentageFeatureSet(gd.singlet.CM, pattern = "^mt-")
VlnPlot(gd.singlet.CM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = T)
gd.singlet.CM <- FindVariableFeatures(gd.singlet.CM, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(gd.singlet.CM), 10)
plot1 <- VariableFeaturePlot(gd.singlet.CM)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
gd.singlet.CM <- ScaleData(gd.singlet.CM)
gd.singlet.CM <- RunPCA(gd.singlet.CM, features = VariableFeatures(object = gd.singlet.CM))
print(gd.singlet.CM[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gd.singlet.CM, dims = 1:2, reduction = "pca")
gd.singlet.CM <- FindNeighbors(gd.singlet.CM, dims = 1:30)
gd.singlet.CM <- FindClusters(gd.singlet.CM)
gd.singlet.CM <- RunUMAP(gd.singlet.CM, dims = 1:30)

DimPlot(gd.singlet.CM, reduction = "umap", label = T)
DimPlot(gd.singlet.CM, reduction = "umap", label = T,         cells = names(gd.singlet.CM$types)[grep("CMO301", gd.singlet.CM$types)])
DimPlot(gd.singlet.CM, reduction = "umap", label = T,         cells = names(gd.singlet.CM$types)[grep("CMO302", gd.singlet.CM$types)])
DimPlot(gd.singlet.CM, reduction = "umap", label = T,         cells = names(gd.singlet.CM$types)[grep("CMO303", gd.singlet.CM$types)])
DimPlot(gd.singlet.CM, reduction = "umap", label = T,         cells = names(gd.singlet.CM$types)[grep("CMO304", gd.singlet.CM$types)])
DimPlot(gd.singlet.CM, reduction = "umap", label = T,         cells = names(gd.singlet.CM$types)[grep("CMO305", gd.singlet.CM$types)])

VlnPlot(gd.singlet.CM, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, log = T)
VlnPlot(gd.singlet.CM, features = c("nCount_RNA"), pt.size = 0.5)
VlnPlot(gd.singlet.CM, features = c("nFeature_RNA"), pt.size = 0.5)

set.seed(123)
clust.col = sample(rainbow(25))
DimPlot(gd.singlet.CM, reduction = "umap", label = T,cols = clust.col)
DimPlot(gd.singlet.CM, group.by = "HTO_classification",reduction = "umap",cols = c("#FF99CC", "#CC6699",      "#FFCC00", "#FF9900",     "#CCFF00","#669900",    "#66FFFF","#33CCFF",     "#CC99FF","#6633FF"))

FeaturePlot(gd.singlet.CM, features = c("Pdgfra","Pecam1","Ptprc", "Mcam"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet.CM, features = c("hto_CMO301", "hto_CMO302", "hto_CMO303", "hto_CMO304", "hto_CMO305"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet.CM, features = c("Lifr", "Il6st", "Plxnb1", "Plxna4", "Plxna2", "Erbb4", "Ankr3", "Notch2", "Pard3", "Ednra", "Dsg2", "Cdh4", "Cdh2"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet.CM, features = c("Nkx2-5", "Sox4", "Cenpa", "Tmsb10", "Zfp706"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet.CM, features = c("Ctla4", "Mb", "Mdh2", "Atp5b", "S100a6", "Erbb2", "Tnni1", "Actc1"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet.CM, features = c("Actn2", "Cacna1c", "Kcnj3", "Scn5a","Myh6", "Tnnt2"), min.cutoff = "q9", cols = c("grey90","red4"))
FeaturePlot(gd.singlet.CM, features = c("Erbb2"), min.cutoff = "q9", cols = c("grey90","red4"))

FeaturePlot(gd.singlet.CM, features = c("Ankrd1"), min.cutoff = "q9", cols = c("grey90","red4"))

dediff_marker_gene_list <- list(c("Ctla4", "Mb", "Mdh2", "Atp5b", "Nkx2-5", "Sox4", "Cenpa", "S100a6", "Erbb2", "Tnni1", "Actc1", "Tmsb10", "Zfp706"))
dediff_marker_TFgene_list <- list(c("Nkx2-5", "Sox4", "Cenpa", "Tmsb10", "Zfp706"))
cycling_marker_gene_list <- list(c("Mki67", "Ccnd3", "Cdk14", "Knl1", "Kif11", "Cdk14", "Pcna"))
ACM_marker_gene_list <- list(c("Actn2", "Cacna1c", "Kcnj3", "Myh6", "Scn5a", "Tnnt2"))

gd.singlet.CM <- AddModuleScore(object = gd.singlet.CM, features = dediff_marker_gene_list, name = "Dedifferentiation_score")
gd.singlet.CM <- AddModuleScore(object = gd.singlet.CM, features = dediff_marker_TFgene_list, name = "Dedifferentiation_TF_score")
gd.singlet.CM <- AddModuleScore(object = gd.singlet.CM, features = cycling_marker_gene_list, name = "Cycling_score")
gd.singlet.CM <- AddModuleScore(object = gd.singlet.CM, features = ACM_marker_gene_list, name = "ACM_score")
FeaturePlot(object = gd.singlet.CM, features = "ACM_score1")
FeaturePlot(object = gd.singlet.CM, features = "Dedifferentiation_score1")
FeaturePlot(object = gd.singlet.CM, features = "Dedifferentiation_TF_score1")
FeaturePlot(object = gd.singlet.CM, features = "Cycling_score1")

table(gd.singlet.CM@meta.data$HTO_classification)


# find markers for every cluster compared to all remaining cells, report only the positive ones
library(tidyr)
library(dplyr)
gd.singlet.CM.markers <- FindAllMarkers(gd.singlet.CM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gd.singlet.CM.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(gd.singlet.CM, features = top10$gene) + NoLegend()


# Visualize co-expression of two features simultaneously
FeaturePlot(gd.singlet.CM, features = c("Actc1", "Atp5b"), blend = TRUE)
FeaturePlot(gd.singlet.CM, features = c("Mb", "Mdh2"), blend = TRUE)
FeaturePlot(gd.singlet.CM, features = c("Tnnt2", "Scn5a"), blend = TRUE)
FeaturePlot(gd.singlet.CM, features = c("Myh6", "Tnnt2"), blend = TRUE)





# list having all information of singlet names and CMO identities:
head(gd.singlet.CM@meta.data[,9])
head(rownames(gd.singlet.CM@meta.data))
# collect singlet cells
SingletsCMOlist <- gd.singlet.CM@meta.data[,9]
names(SingletsCMOlist) <- rownames(gd.singlet.CM@meta.data)
length(SingletsCMOlist)
head(SingletsCMOlist)
table(SingletsCMOlist)

# save each time point as 1 separate .rds
heartSinglets <- names(SingletsCMOlist)[SingletsCMOlist%in%c("CMO305")]
length(heartSinglets)

dim(gd.umis)
prdata_MI_CM <-gd.umis[,heartSinglets]
dim(prdata_MI_CM)
saveRDS(prdata_MI_CM, file = "MI_CM_expmatrix.rds")
