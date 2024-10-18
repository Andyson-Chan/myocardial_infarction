
library(DOSE)
library(enrichplot)
library(ReactomePA)

plotmap(sc1000_LAD_Cryo_Sham_neural, um=T)
table(sc1000_LAD_Cryo_Sham_neural@cpart)
# compare difference between Quiescent and Galectin SwC
A <- names(sc1000_LAD_Cryo_Sham_neural@cpart)[sc1000_LAD_Cryo_Sham_neural@cpart %in% c(1)]
B <- names(sc1000_LAD_Cryo_Sham_neural@cpart)[sc1000_LAD_Cryo_Sham_neural@cpart %in% c(2)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_neural,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Quiescent SwC",Bname="Galectin SwC",show_names=TRUE,padj=TRUE)
plotexpmap(sc1000_LAD_Cryo_Sham_neural, logsc=T, um=T, "Tubb2a", cex=2)

SwC_DEGs_FC<- x$res[order(-x$res$log2FoldChange),]
entrez <- bitr(rownames(SwC_DEGs_FC), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
SwC_DEG_entrez <- SwC_DEGs_FC[entrez$SYMBOL,"log2FoldChange"]
names(SwC_DEG_entrez)<-entrez$ENTREZID
head(SwC_DEG_entrez)

geneList <- SwC_DEG_entrez
length(geneList)



# Gene-Concept Network
MP_df <- data.frame(Entrez=names(geneList), log2FC=geneList)
MP_df <- MP_df[abs(MP_df$log2FC) > 0.3,]
MP_df$group <- "Galectin SwC"
MP_df$group[MP_df$log2FC < -0.3] <- "Quiescent SwC"

# create a list summerizing for different samples and their enriched genes (needed for some plots)
list<-list(MP_df$Entrez[MP_df$group%in%"Galectin SwC"],MP_df$Entrez[MP_df$group%in%"Quiescent SwC"])
names(list)<-c("Galectin SwC","Quiescent SwC")

# make and plot gene concept network (GCK) and pathway dotplots
# database selection can be from One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" (only some works in mouse)

# Database: enrichPathway
GCN <- compareCluster(geneCluster = list, fun = enrichPathway, organism = "mouse") 
GCN@readable <- FALSE
GCN <- setReadable(GCN, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(GCN ,cex_label_gene = 1,cex_label_category= 1)
