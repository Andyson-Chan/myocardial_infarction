
library(DOSE)
library(enrichplot)
library(ReactomePA)


# compare difference between infiltrating and resident macrophages
A <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart %in% c(8,7,6,13)]
B <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart %in% c(2,1,20,22)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_Immune1,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="CCR2+ Macrophages",Bname="Resident Macrophages",show_names=TRUE,padj=TRUE)

# prepare geneList (contains geneIDs and log2FC values)
macro_DEGs_FC<- x$res[order(-x$res$log2FoldChange),]
entrez <- bitr(rownames(macro_DEGs_FC), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
macro_DEGs_entrez <- macro_DEGs_FC[entrez$SYMBOL,"log2FoldChange"]
names(macro_DEGs_entrez)<-entrez$ENTREZID
head(macro_DEGs_entrez)

geneList <- macro_DEGs_entrez


# select genes from either populations

# collect genes enriched in resident macrophages (log2FC > 0)
de <- names(macro_DEGs_entrez)[geneList > 0.5]
# collect genes enriched in CCR2+ macrophages (log2FC < 0)
de <- names(macro_DEGs_entrez)[geneList < -0.5]


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



# Gene-Concept Network
MP_df <- data.frame(Entrez=names(geneList), log2FC=geneList)
MP_df <- MP_df[abs(MP_df$log2FC) > 1,]
MP_df$group <- "Resident MP"
MP_df$group[MP_df$log2FC < 1] <- "Circulating MP"

# create a list summerizing for different samples and their enriched genes (needed for some plots)
list<-list(MP_df$Entrez[MP_df$group%in%"Resident MP"],MP_df$Entrez[MP_df$group%in%"Circulating MP"])
names(list)<-c("Resident MP","Circulating MP")

# make and plot gene concept network (GCK) and pathway dotplots
# database selection can be from One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" (only some works in mouse)

# Database: enrichPathway
GCN <- compareCluster(geneCluster = list, fun = enrichPathway, organism = "mouse") 
GCN@readable <- FALSE
GCN <- setReadable(GCN, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(GCN ,cex_label_gene = 1,cex_label_category= 1)

formula_res <- compareCluster(Entrez~group, data=MP_df, fun = enrichPathway , organism = "mouse")
dotplot(formula_res)


# Heatmap-like functional classification
heatplot(edox, foldChange=geneList, showCategory=5)


# Enrichment Map
edo <- pairwise_termsim(edo)
emapplot(edo)



# prepare GSEA object

# from GO terms
gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")

gse$Description[1:50]

# from Reactome
gse <- gsePathway(geneList=geneList, 
                  minGSSize = 15, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = F, 
                  organism="mouse", 
                  pAdjustMethod = "BH")

gse$Description[1:50]



# ridgeline plot for expression distribution of GSEA result

ridgeplot(gse)+ labs(x = "enrichment distribution")

# running score and preranked list of GSEA result (pick the relevent pathways)
require(enrichplot)
gseaplot2(gse, geneSetID = c(2,3,4,6,9,17,27,32,36,40,41,42))

# search for the corresponding geneSetIDs
gse@result[1:3,c("ID","Description")]

geneSetID <- c(
  rownames(gse@result)[gse@result$Description%in%"Complement cascade"],
  rownames(gse@result)[gse@result$Description%in%"Vesicle-mediated transport"],
  rownames(gse@result)[gse@result$Description%in%"Cellular response to hypoxia"],
  rownames(gse@result)[gse@result$Description%in%"Neutrophil degranulation"]
)
gseaplot2(gse, geneSetID = geneSetID)
gseaplot2(gse, geneSetID = c("Innate Immune System", "Translation"))






# list of Entrez gene IDs that contributed to the enrichment score
genelist <- fgseaRes[fgseaRes$pathway%in%"Complement cascade",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist


pdf("/Users/andychan/Desktop/Deep_sequencing/adult_all_datasets_merged/panels/immune cells/umap_Immune_Complementcascade.pdf", width=7.5, height=5.5) #defaults are 7 and 7
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, genelist$SYMBOL, n = "Complement cascade", um=T, cex=0.3, logsc=T)
dev.off()

pdf("/Users/andychan/Desktop/Deep_sequencing/adult_all_datasets_merged/panels/immune cells/umap_Immune_Proinf_cytokines.pdf", width=7.5, height=5.5) #defaults are 7 and 7
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, cytokines_proinf, n = "Proinflammatory cytokine score", um=T, cex=0.3, logsc=T)
dev.off()


genelist <- fgseaRes[fgseaRes$pathway%in%"Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

pdf("/Users/andychan/Desktop/Deep_sequencing/adult_all_datasets_merged/panels/immune cells/umap_Immune_IGF_transport.pdf", width=7.5, height=5.5) #defaults are 7 and 7
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, genelist$SYMBOL, n = "IGF transport and uptake by IGFBPs", um=T, cex=0.3, logsc=T)
dev.off()

genelist <- fgseaRes[fgseaRes$pathway%in%"Cellular response to hypoxia",leadingEdge]
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist

pdf("/Users/andychan/Desktop/Deep_sequencing/adult_all_datasets_merged/panels/immune cells/umap_Immune_hypoxia.pdf", width=7.5, height=5.5) #defaults are 7 and 7
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, genelist$SYMBOL, n = "Cellular response to hypoxia", um=T, cex=0.3, logsc=T)
dev.off()

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