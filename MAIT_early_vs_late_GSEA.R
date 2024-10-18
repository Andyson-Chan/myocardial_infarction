# general MAIT-like cell population dynamics (% of all immune cells)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]

MAIT_cells <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(11,18)]
v1 <- length(intersect(sample1, MAIT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v3 <- length(intersect(sample3, MAIT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v7 <- length(intersect(sample7, MAIT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v28 <- length(intersect(sample28, MAIT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v56 <- length(intersect(sample56, MAIT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
vSh <- length(intersect(sampleSh, MAIT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_Sh|NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
a<-c(vSh,v1,v3,v7,v28,v56)
names(a)<-timepoints
barplot(a, xlab = "Days post-MI", ylab = "percentage of MAIT-like cells",main = "% of MAIT-like cells in immune compartment", col = c("#CC33FF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000CC"))




library(DOSE)
library(enrichplot)
library(ReactomePA)

plotmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, um=T)
table(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)
# compare difference between early and late MAIT-like cells
A <- intersect(c(sample1, sample3, sample7), MAIT_cells)
B <- intersect(c(sample28, sample56), MAIT_cells)
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Early MAIT-like cells",Bname="Late MAIT-like cells",show_names=TRUE,padj=TRUE)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, logsc=T, um=T, "Atp5k")

MAIT_DEGs_FC<- x$res[order(-x$res$log2FoldChange),]
entrez <- bitr(rownames(MAIT_DEGs_FC), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
MAIT_DEG_entrez <- MAIT_DEGs_FC[entrez$SYMBOL,"log2FoldChange"]
names(MAIT_DEG_entrez)<-entrez$ENTREZID
head(MAIT_DEG_entrez)

geneList <- MAIT_DEG_entrez
length(geneList)

# select genes from either populations

# collect genes enriched in early MAIT-like cells (log2FC > 0)
de <- names(MAIT_DEG_entrez)[geneList > 0]
# collect genes enriched in late MAIT-like cells (log2FC < 0)
de <- names(MAIT_DEG_entrez)[geneList < 0]


# choose the database you want for the enrichment measurements

# Reactome
edo <- enrichPathway(gene=de, organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable=TRUE)
edo <- enrichPathway(gene=de, organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 1, readable=TRUE)

# plot the selected database:
barplot(edo, showCategory=10) 

# Heatmap-like functional classification
edox <- setReadable(edo, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
heatplot(edox, foldChange=geneList, showCategory=4)


# Gene-Concept Network
MP_df <- data.frame(Entrez=names(geneList), log2FC=geneList)
MP_df <- MP_df[abs(MP_df$log2FC) > 0.3,]
MP_df$group <- "Late MAIT-like cells"
MP_df$group[MP_df$log2FC < -0.3] <- "Early MAIT-like cells"

# create a list summerizing for different samples and their enriched genes (needed for some plots)
list<-list(MP_df$Entrez[MP_df$group%in%"Late MAIT-like cells"],MP_df$Entrez[MP_df$group%in%"Early MAIT-like cells"])
names(list)<-c("Late MAIT-like cells","Early MAIT-like cells")

# make and plot gene concept network (GCK) and pathway dotplots
# database selection can be from One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" (only some works in mouse)

# Database: enrichPathway
GCN <- compareCluster(geneCluster = list, fun = enrichPathway, organism = "mouse") 
GCN@readable <- FALSE
GCN <- setReadable(GCN, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(GCN ,cex_label_gene = 1,cex_label_category= 1)


plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, logsc=T, um=T,cex=1, "Bcl2", cells = intersect(c(sample1, sample3, sample7), MAIT_cells))
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, logsc=T, um=T,cex=1, "Bcl2", cells = intersect(c(sample28, sample56), MAIT_cells))


formula_res <- compareCluster(Entrez~group, data=MP_df, fun = enrichPathway , organism = "mouse")
dotplot(formula_res)




# Enrichment Map
edo <- pairwise_termsim(edo)
emapplot(edo)


