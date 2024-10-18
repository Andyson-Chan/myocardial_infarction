
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)

plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, cex=0.2)
# compare difference between cryo and LAD FB
# all myoFB
myoFb_cells <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)[sc1000_LAD_Cryo_Sham_Fibro@cpart %in% c(5,7,17,3,8,15,10,4)]
typesC <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata), 7,8)
sampleCryo <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("MI_",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sampleLAD <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("LA_",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
myoFB_cryo <- intersect(myoFb_cells, sampleCryo)
myoFB_LAD <- intersect(myoFb_cells, sampleLAD)
length(myoFB_cryo); length(myoFB_LAD)


x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_Fibro,n=c(myoFB_cryo,myoFB_LAD)), A=myoFB_cryo, B=myoFB_LAD )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Cryo myoFB",Bname="LAD myoFB",show_names=TRUE,padj=TRUE)

# prepare geneList (contains geneIDs and log2FC values)
FB_DEGs_FC<- x$res[order(-x$res$log2FoldChange),]
entrez <- bitr(rownames(FB_DEGs_FC), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
FB_DEGs_entrez <- FB_DEGs_FC[entrez$SYMBOL,"log2FoldChange"]
names(FB_DEGs_entrez)<-entrez$ENTREZID
head(FB_DEGs_entrez)

geneList <- FB_DEGs_entrez


# select genes from either populations

# collect genes enriched in resident macrophages (log2FC > 0)
de <- names(FB_DEGs_entrez)[geneList > 0.5]
# collect genes enriched in CCR2+ macrophages (log2FC < 0)
de <- names(FB_DEGs_entrez)[geneList < -0.5]


# Reactome database for the enrichment measurements

# Reactome
edo <- enrichPathway(gene=de, organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable=TRUE)

# plot the selected database:
barplot(edo, showCategory=20) 



# Gene-Concept Network
FB_df <- data.frame(Entrez=names(geneList), log2FC=geneList)
FB_df <- FB_df[abs(FB_df$log2FC) > 0.43,]
FB_df$group <- "LAD myoFB"
FB_df$group[FB_df$log2FC < -0.43] <- "Cryo myoFB"

# create a list summerizing for different samples and their enriched genes (needed for some plots)
list<-list(FB_df$Entrez[FB_df$group%in%"LAD myoFB"],FB_df$Entrez[FB_df$group%in%"Cryo myoFB"])
names(list)<-c("LAD myoFB","Cryo myoFB")

# make and plot gene concept network (GCK) and pathway dotplots
# database selection can be from One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" (only some works in mouse)

# Database: enrichPathway
GCN <- compareCluster(geneCluster = list, fun = enrichPathway, organism = "mouse") 
GCN@readable <- F
GCN <- setReadable(GCN, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(GCN ,cex_label_gene = 1,cex_label_category= 1)

formula_res <- compareCluster(Entrez~group, data=FB_df, fun = enrichPathway , organism = "mouse")
dotplot(formula_res)


# Heatmap-like functional classification
heatplot(GCN, foldChange=geneList, showCategory=5)


# Enrichment Map
edo <- pairwise_termsim(edo)
emapplot(edo)
