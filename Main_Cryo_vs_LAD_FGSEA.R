# compare DEG between LAD and Cryoablated cells
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 7,8)
head(types)
unique(types)
types[types%in%"LA"] <- "LAD"
types[types%in%"MI"] <- "Cryoablation"
types[types%in%"Sh"] <- "Sham"
names(types) <- colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)
A <- names(types[types%in%"LAD"]) # LAD
B <- names(types[types%in%"Cryoablation"]) # Cryo
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_harmony,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="LAD",Bname="Cryoablation",show_names=TRUE,padj=TRUE)




library(fgsea)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)


# arrange fold changes from small to large values (-)
DEGs<-x$res[order(-x$res$log2FoldChange),]
#DEGs<-DEGs[cl1vs8$padj<0.05,]
entrez <- bitr(rownames(DEGs), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
#cl1vs8_entrez<-data.frame(entrez$ENTREZID,cl1vs8[entrez$SYMBOL,"log2FoldChange"])

DEGs_entrez <- DEGs[entrez$SYMBOL,"log2FoldChange"]
names(DEGs_entrez) <- entrez$ENTREZID
barplot(sort(DEGs_entrez, decreasing = T))

data(examplePathways)
library(reactome.db)
REpathways <- reactomePathways(names(DEGs_entrez))


# there are 2 databases available for choice: REpathways, and examplePathways

#p-value less reliable
fgseaRes <- fgsea(pathways = REpathways, 
                  stats    = DEGs_entrez,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes[order(pval), ])

#p-value more accurate
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = DEGs_entrez,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes[order(pval), ])

#gene_list<-entrez[entrez$SYMBOL,fgseaRes$leadingEdge]

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
dev.off()
plotGseaTable(examplePathways[topPathways], DEGs_entrez, fgseaRes, 
              gseaParam=0.5)
plotGseaTable(REpathways[topPathways], DEGs_entrez, fgseaRes, 
              gseaParam=0.5)

# plot a particular pathway enlarged:
plotEnrichment(REpathways[["Cellular response to hypoxia"]],
               DEGs_entrez) + labs(title="Cellular response to hypoxia") 
plotEnrichment(REpathways[["Complement cascade"]],
               DEGs_entrez) + labs(title="Complement cascade") 


dev.off()
#collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
#                                      examplePathways, DEGs_entrez)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][pval < 0.05], 
                                      examplePathways, DEGs_entrez)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(examplePathways[mainPathways], DEGs_entrez, fgseaRes, 
              gseaParam = 0.5)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][pval < 0.05], 
                                      REpathways, DEGs_entrez)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(REpathways[mainPathways], DEGs_entrez, fgseaRes, 
              gseaParam=0.5)


# list of Entrez gene IDs that contributed to the enrichment score
genelist <- fgseaRes[fgseaRes$pathway%in%"Muscle contraction",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

plotexpmap(sc1000_LAD_Cryo_Sham_harmony, genelist$SYMBOL, n = "Muscle contraction in LAD", um=F, cex=0.2, logsc=T, cells = names(types[types%in%"LAD"]))
plotexpmap(sc1000_LAD_Cryo_Sham_harmony, genelist$SYMBOL, n = "Muscle contraction in Cryoablation", um=F, cex=0.2, logsc=T, cells = names(types[types%in%"Cryoablation"]))

genelist <- fgseaRes[fgseaRes$pathway%in%"mRNA Splicing",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

plotexpmap(sc1000_LAD_Cryo_Sham_harmony, genelist$SYMBOL, n = "mRNA Splicing in LAD", um=F, cex=0.2, logsc=T, cells = names(types[types%in%"LAD"]))
plotexpmap(sc1000_LAD_Cryo_Sham_harmony, genelist$SYMBOL, n = "mRNA Splicing in Cryoablation", um=F, cex=0.2, logsc=T, cells = names(types[types%in%"Cryoablation"]))

















sampleLAD <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("LA_",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex=0.5, logsc=T, cells=sampleLAD, genelist$SYMBOL, n = "Cellular response to hypoxia, LAD")
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex=0.5, logsc=T, cells=sampleLAD, "Hif1a", n = "Hif1a expression in LAD")
sampleCryo <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("MI_",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex=0.5, logsc=T, cells=sampleCryo, genelist$SYMBOL, n = "Cellular response to hypoxia, Cryoablation")
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex=0.5, logsc=T, cells=sampleCryo, "Hif1a", n = "Hif1a expression in Cryoablation")


genelist <- fgseaRes[fgseaRes$pathway%in%"Clathrin-mediated endocytosis",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, genelist$SYMBOL, n = "Clathrin-mediated endocytosis", um=T, cex=0.2, logsc=T)


genelist <- fgseaRes[fgseaRes$pathway%in%"5992210_Chondroitin_sulfate_dermatan_sulfate_metabolism",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

genelist <- fgseaRes[fgseaRes$pathway%in%"5991838_Netrin-1_signaling",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

genelist <- fgseaRes[fgseaRes$pathway%in%"5991757_RHO_GTPases_Activate_Formins",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

genelist <- fgseaRes[fgseaRes$pathway%in%"5991956_Cell_junction_organization",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist


