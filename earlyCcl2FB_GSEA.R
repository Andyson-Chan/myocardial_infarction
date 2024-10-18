
library(DOSE)
library(enrichplot)
library(ReactomePA)

plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, cex=0.2)
table(sc1000_LAD_Cryo_Sham_Fibro@cpart)
# check DEGs of early myoFB (cluster 5)

A <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)[sc1000_LAD_Cryo_Sham_Fibro@cpart %in% c(1:4, 6:19)]
B <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)[sc1000_LAD_Cryo_Sham_Fibro@cpart %in% c(5)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_Fibro,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="FB/myoFB",Bname="early myoFB",show_names=TRUE,padj=TRUE)

plotexpmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, logsc=T, "Mgp")
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, logsc=T, "Il2ra")
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, logsc=T, "Ftl1")
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, logsc=T, "Sdhb")
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, logsc=T, "Surf1")

# Collect DEGs of tFB cells
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Fibro,c(5),pvalue=0.03)
head(dg$dg)
FB_DEGs_FC<- dg$dg[order(-dg$dg$fc),]
head(FB_DEGs_FC)
FB_DEGs_FC <- FB_DEGs_FC[,"fc"]
names(FB_DEGs_FC) <- rownames(dg$dg[order(-dg$dg$fc),])
FB_DEGs_FC <- log2(FB_DEGs_FC)

entrez <- bitr(names(FB_DEGs_FC), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
FB_DEGs_entrez <- FB_DEGs_FC[entrez$SYMBOL]
names(FB_DEGs_entrez)<-entrez$ENTREZID
head(FB_DEGs_entrez)

geneList <- FB_DEGs_entrez


# select genes from either populations

# collect genes enriched in tFB (log2FC > 0)
de <- names(FB_DEGs_entrez)[geneList > 0]
# collect genes downregulated in tFB (log2FC < 0)
de <- names(FB_DEGs_entrez)[geneList < -0.5]


# choose the database you want for the enrichment measurements

# Reactome
edo <- enrichPathway(gene=de, organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable=TRUE)
# Wikipathway
edo <- enrichWP(de, organism = 'Mus musculus')
# Go Terms (3 databases)
edo <- enrichGO(geneList, universe=uni[,2], OrgDb = "org.Mm.eg.db", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
edo <- enrichGO(geneList, universe=uni[,2], OrgDb = "org.Mm.eg.db", ont = "MF",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
edo <- enrichGO(geneList, universe=uni[,2], OrgDb = "org.Mm.eg.db", ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
# KEGG
edo <- enrichKEGG(geneList, organism = 'mmu', keyType = "kegg", pvalueCutoff = 0.05, universe=uni[,2])


# plot the selected database:
barplot(edo, showCategory=20) 


# Heatmap-like functional classification
edox <- setReadable(edo, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
heatplot(edox, foldChange=geneList, showCategory=2)



# Gene-Concept Network
FB_df <- data.frame(Entrez=names(geneList), log2FC=geneList)
FB_df <- FB_df[abs(FB_df$log2FC) > 1,]
FB_df$group <- "early myoFB"
FB_df$group[FB_df$log2FC < 1] <- "FB/myoFB"

# create a list summerizing for different samples and their enriched genes (needed for some plots)
list<-list(FB_df$Entrez[FB_df$group%in%"LAD myoFB"],FB_df$Entrez[FB_df$group%in%"Cryo myoFB"])
names(list)<-c("tFB","FB/myoFB")

# make and plot gene concept network (GCK) and pathway dotplots
# database selection can be from One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" (only some works in mouse)

# Database: enrichPathway
GCN <- compareCluster(geneCluster = list, fun = enrichPathway, organism = "mouse") 
GCN@readable <- F
GCN <- setReadable(GCN, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(GCN ,cex_label_gene = 1,cex_label_category= 1)

formula_res <- compareCluster(Entrez~group, data=FB_df, fun = enrichPathway , organism = "mouse")
dotplot(formula_res)


