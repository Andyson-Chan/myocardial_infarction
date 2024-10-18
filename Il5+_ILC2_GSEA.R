
library(DOSE)
library(enrichplot)
library(ReactomePA)

plotmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, um=T)
table(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)
# compare difference between IL5+/- ILC2
A <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart %in% c(1)]
B <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart %in% c(17)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Il5- ILC2",Bname="Il5+ ILC2",show_names=TRUE,padj=TRUE)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, logsc=T, um=T, "Cd3e")

ILC2_DEGs_FC<- x$res[order(-x$res$log2FoldChange),]
entrez <- bitr(rownames(ILC2_DEGs_FC), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
ILC2_DEG_entrez <- ILC2_DEGs_FC[entrez$SYMBOL,"log2FoldChange"]
names(ILC2_DEG_entrez)<-entrez$ENTREZID
head(ILC2_DEG_entrez)

geneList <- ILC2_DEG_entrez
length(geneList)

# select genes from either populations

# collect genes enriched in Il5+ ILC2 (log2FC > 0)
de <- names(ILC2_DEG_entrez)[geneList > 0]
# or collect genes enriched in IL5- ILC2 (log2FC < 0)
de <- names(ILC2_DEG_entrez)[geneList < 0]


# choose the database you want for the enrichment measurements

# Reactome
edo <- enrichPathway(gene=de, organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable=TRUE)
edo <- enrichPathway(gene=de, organism="mouse", pAdjustMethod = "BH", pvalueCutoff = 1, readable=TRUE)
# Wikipathway
edo <- enrichWP(de, organism = 'Mus musculus')
# Go Terms (3 databases)
edo <- enrichGO(geneList, universe=uni[,2], OrgDb = "org.Mm.eg.db", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
edo <- enrichGO(geneList, universe=uni[,2], OrgDb = "org.Mm.eg.db", ont = "MF",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
edo <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
# KEGG
edo <- enrichKEGG(geneList, organism = 'mmu', keyType = "kegg", pvalueCutoff = 0.05, universe=uni[,2])


# plot the selected database:
barplot(edo, showCategory=10) 

# Heatmap-like functional classification
edox <- setReadable(edo, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
heatplot(edox, foldChange=geneList, showCategory=4)


# Gene-Concept Network
MP_df <- data.frame(Entrez=names(geneList), log2FC=geneList)
MP_df <- MP_df[abs(MP_df$log2FC) > 0.3,]
MP_df$group <- "Il5+ ILC2"
MP_df$group[MP_df$log2FC < -0.3] <- "Il5- ILC2"

# create a list summerizing for different samples and their enriched genes (needed for some plots)
list<-list(MP_df$Entrez[MP_df$group%in%"Il5+ ILC2"],MP_df$Entrez[MP_df$group%in%"Il5- ILC2"])
names(list)<-c("Il5+ ILC2","Il5- ILC2")

# make and plot gene concept network (GCK) and pathway dotplots
# database selection can be from One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" (only some works in mouse)

# Database: enrichPathway
GCN <- compareCluster(geneCluster = list, fun = enrichPathway, organism = "mouse") 
GCN@readable <- FALSE
GCN <- setReadable(GCN, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(GCN ,cex_label_gene = 1,cex_label_category= 1)

formula_res <- compareCluster(Entrez~group, data=MP_df, fun = enrichPathway , organism = "mouse")
dotplot(formula_res)




# Enrichment Map
edo <- pairwise_termsim(edo)
emapplot(edo)

