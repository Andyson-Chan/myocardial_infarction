# subcluster fibroblasts 

Fibroblasts <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)[sc1000_LAD_Cryo_Sham_harmony@cpart%in%c(17,21,2,13,26,3,1,30)]


sc1000_LAD_Cryo_Sham_Fibro <-SCseq(sc1000_LAD_Cryo_Sham_harmony@expdata[,Fibroblasts])
sc1000_LAD_Cryo_Sham_Fibro<-filterdata(sc1000_LAD_Cryo_Sham_Fibro,mintotal=1000, FGenes=rownames(sc1000_LAD_Cryo_Sham_Fibro@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(sc1000_LAD_Cryo_Sham_Fibro@expdata))])
expData  <- getExpData(sc1000_LAD_Cryo_Sham_Fibro)
# with batch removal (harmony)
batches <- Fibroblasts
names(batches) <- Fibroblasts
batch1 <- batches[grep("D01|D03|D07|D28|MI_D56",batches)]
batch2 <- batches[grep("LA_D56|Sh_D56",batches)]
batch1 <- replace(batch1,,"b1")
batch2 <- replace(batch2,,"b2")
batchesFB <- c(batch1, batch2)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesFB, bmethod = "harmony")
# no batch effect removal
# res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham_Fibro@expdata)[sc1000_LAD_Cryo_Sham_Fibro@genes,colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)]
noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=10)
sc1000_LAD_Cryo_Sham_Fibro <- updateSC(sc1000_LAD_Cryo_Sham_Fibro,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham_Fibro <- comptsne(sc1000_LAD_Cryo_Sham_Fibro)
sc1000_LAD_Cryo_Sham_Fibro <- compumap(sc1000_LAD_Cryo_Sham_Fibro, spread = 1, min_dist = 0.5, n_neighbors = 30)
plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, cex = 0.2)
plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=F, cex = 0.2)
table(sc1000_LAD_Cryo_Sham_Fibro@cpart)



batches <- Fibroblasts
names(batches) <- Fibroblasts
batch1 <- batches[grep("D01|D03|D07|D28|MI_D56",batches)]
batch2 <- batches[grep("LA_D56|Sh_D56",batches)]
batch1 <- replace(batch1,,"b1")
batch2 <- replace(batch2,,"b2")
batchesFB <- c(batch1, batch2)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesFB, bmethod = "harmony")
# no batch effect removal
# res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham_Fibro@expdata)[sc1000_LAD_Cryo_Sham_Fibro@genes,colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)]
noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=10)
sc1000_LAD_Cryo_Sham_Fibro1 <- updateSC(sc1000_LAD_Cryo_Sham_Fibro,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham_Fibro1 <- comptsne(sc1000_LAD_Cryo_Sham_Fibro1)
sc1000_LAD_Cryo_Sham_Fibro1 <- compumap(sc1000_LAD_Cryo_Sham_Fibro1, spread = 1, min_dist = 0.5, n_neighbors = 30)
plotmap(sc1000_LAD_Cryo_Sham_Fibro1, um=T, cex = 0.2)
plotmap(sc1000_LAD_Cryo_Sham_Fibro1, um=F, cex = 0.2)
table(sc1000_LAD_Cryo_Sham_Fibro1@cpart)










# change cluster colors
library(ggplot2)
library(paletteer)
sc1000_LAD_Cryo_Sham_Fibro_fcol_original <- sc1000_LAD_Cryo_Sham_Fibro@fcol

Fibro_fcol <- sc1000_LAD_Cryo_Sham_Fibro@fcol
Fibro_fcol[c(1,2,6,19,12,13)] <- as.character(paletteer_c("grDevices::Spectral", 30))[c(25:30)]
Fibro_fcol[c(5,14,9)] <- as.character(paletteer_c("grDevices::Spectral", 30))[c(10,17,20)]
Fibro_fcol[c(15,10,8,17,7,4)] <- as.character(paletteer_c("grDevices::Spectral", 30))[c(1:6)]
Fibro_fcol[c(16,11,18)] <- as.character(paletteer_c("grDevices::RdGy", 3))
sc1000_LAD_Cryo_Sham_Fibro@fcol <- Fibro_fcol

plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, cex=0.5)


# plot umap with cluster names annotated

cluster_ann <- sc1000_LAD_Cryo_Sham_Fibro@cpart
cluster_ann[cluster_ann%in%1] <- "FB_quiescent"
cluster_ann[cluster_ann%in%2] <- "FB_quiescent"
cluster_ann[cluster_ann%in%3] <- "FBmyo_IFN" 
cluster_ann[cluster_ann%in%4] <- "FB_Postn_Thbs4" 
cluster_ann[cluster_ann%in%5] <- "FB_Ccl2" 
cluster_ann[cluster_ann%in%6] <- "FB_Duox1"
cluster_ann[cluster_ann%in%7] <- "FBmyo" 
cluster_ann[cluster_ann%in%8] <- "FBmyo_IFN"
cluster_ann[cluster_ann%in%9] <- "FB_transition_Cd9"
cluster_ann[cluster_ann%in%10] <- "Fbmyo_Dkk2"
cluster_ann[cluster_ann%in%11] <- "FB_rMacro" 
cluster_ann[cluster_ann%in%12] <- "FB_Fgl2"
cluster_ann[cluster_ann%in%13] <- "FB_Cxcl14" 
cluster_ann[cluster_ann%in%14] <- "FB_transition_Il1rapl1"
cluster_ann[cluster_ann%in%15] <- "Fbmyo_cycling"
cluster_ann[cluster_ann%in%16] <- "FB_Macro_Saa3" 
cluster_ann[cluster_ann%in%17] <- "Fbmyo_Hp" 
cluster_ann[cluster_ann%in%18] <- "FB_Neutrophil" 
cluster_ann[cluster_ann%in%19] <- "FB_WntX" 

cluster_ann_Fibro <- cluster_ann

plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, cluster_ann_Fibro, um=T, cex=0.2, samples_col = 
                 as.character(paletteer_c("grDevices::Spectral", 19))
)






# visualize different conditions
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
typesC <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata), 7,8)
head(typesC)
unique(typesC)
typesC[typesC%in%"LA"] <- "LAD"
typesC[typesC%in%"MI"] <- "Cryoablation"
typesC[typesC%in%"Sh"] <- "Sham"
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesC, um=T, cex = 0.2, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesC, um=T, cex = 0.2, subset = "LAD")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesC, um=T, cex = 0.2, subset = "Cryoablation")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesC, um=T, subset = "Sham", samples_col  = rep("#0A1722", length(typesC[typesC%in%"Sham"])))

# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
typesT <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata), 7,12)
# merge all Sham time points into 1 type
typesT<- gsub("Sh_D01", "Sham", typesT)
typesT<- gsub("Sh_D03", "Sham", typesT)
# typesT<- gsub("Sh_D07", "Sham", typesT) contaminated with LAD cells
typesT<- gsub("Sh_D28", "Sham", typesT)
typesT<- gsub("Sh_D56", "Sham", typesT)

typesT<- gsub("LA_D01", "D01", typesT)
typesT<- gsub("LA_D03", "D03", typesT)
typesT<- gsub("LA_D07", "D07", typesT)
typesT<- gsub("LA_D28", "D28", typesT)
typesT<- gsub("LA_D56", "D56", typesT)

typesT<- gsub("MI_D01", "D01", typesT)
typesT<- gsub("MI_D03", "D03", typesT)
typesT<- gsub("MI_D07", "D07", typesT)
typesT<- gsub("MI_D28", "D28", typesT)
typesT<- gsub("MI_D56", "D56", typesT)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=T, cex = 0.2, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=F, subset = "D01")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=F, subset = "D03")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=T, subset = "D07", leg=F, samples_col  = rep("#0E7535", length(typesT[typesT%in%"D07"])))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=T, subset = "D28")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=T, subset = "D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=T, cex=0.2, subset = c("Sham", "D01","D03"),samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Fibro, typesT, um=T, cex=0.2, subset = c("Sham", "D07","D28","D56"),samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))

plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Pdgfra",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Acta2",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Postn",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Thbs4",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Il6",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Eln",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Wfdc18",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Saa3",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Hp",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Ccl2",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Cxcl5",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Il11",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Cxcl3",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Mki67",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Ptprc",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"S100a9",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Itgam",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Pdgfc",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Col11a1",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Dkk2",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Cd109",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Isg15",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Cxcl10",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Tnc",logsc=T,fr=F, cex = 0.5, um=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Mfap5",logsc=T,fr=F, cex = 0.5, um=T)

A <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)[sc1000_LAD_Cryo_Sham_Fibro@cpart %in% c(4)]
B <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)[sc1000_LAD_Cryo_Sham_Fibro@cpart %in% c(10)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_Fibro,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="D7 myoFB LAD+Cryo",Bname="D7 myoFB Cryo-specific",show_names=TRUE,padj=TRUE)

plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Cdh5",logsc=T,fr=F, um=T)

dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Fibro,c(9),pvalue=0.01)

# select and assign new cluster, medoids and color for Wif1+ Erbb4+ FB:

require(Seurat)
library(patchwork)
library(ggplot2)

plot<-ggplot(sc1000_LAD_Cryo_Sham_Fibro@umap, aes(V1, V2)) +
  geom_point()
select.cells <- CellSelector(plot = plot)


cpart <- sc1000_LAD_Cryo_Sham_Fibro@cpart
cpart[select.cells] <- max(sc1000_LAD_Cryo_Sham_Fibro@cpart) + 1
unique(cpart)
head(cpart)
head(cpart[select.cells])
m <- compmedoids(sc1000_LAD_Cryo_Sham_Fibro,cpart)
fcol <- sc1000_LAD_Cryo_Sham_Fibro@fcol
fcol[length(fcol)+1] <- "#1F5591"
sc1000_LAD_Cryo_Sham_Fibro@cpart <- cpart
sc1000_LAD_Cryo_Sham_Fibro@fcol <- fcol
sc1000_LAD_Cryo_Sham_Fibro@medoids <- m
plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, cex = 0.5)
table(sc1000_LAD_Cryo_Sham_Fibro@cpart)










# dotplots of selected genes
# pool clusters into defined sub-populations
plotmap(sc1000_LAD_Cryo_Sham_Fibro, um=T, cex=0.5)
FB_populations <- cluster_ann_Fibro
table(FB_populations)

# matching cell type names with colors
sc1000_LAD_Cryo_Sham_Fibro@fcol
head(sc1000_LAD_Cryo_Sham_Fibro@cpart)

Fibro_fcol <- sc1000_LAD_Cryo_Sham_Fibro@fcol[sc1000_LAD_Cryo_Sham_Fibro@cpart]
names(Fibro_fcol) <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)
head(Fibro_fcol)

head(FB_populations)
identical(names(FB_populations), names(Fibro_fcol))

Fibro_clusters_fcol <- data.frame(clusters=FB_populations, colors=Fibro_fcol)

saveRDS(Fibro_clusters_fcol, "df_cluster_colors_FB_new.rds")

# interested genes
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Fibro,c(2,3),pvalue=0.01)

genes <- c("Cxcl14","Ackr2", "Duox1", "Cmah", "Fgl2", "Nrxn1","Wif1", "Dkk3", "Ccl2", "Cxcl3", 
           "Acta2", "Mki67", "Pcna", "Saa3","Hp","Isg15", "Cd109", "Col11a1", "Postn", "Thbs4", 
           "Cd9", "Cd302", "Il1rapl1"
)

fractDotPlot(sc1000_LAD_Cryo_Sham_Fibro, genes, samples = FB_populations, subset=c(
  "FB_Cxcl14", "FB_Duox1", "FB_Fgl2", "FB_WntX",
  "FB_Ccl2", "Fbmyo_cycling", "Fbmyo_Dkk2", "Fbmyo_Hp", "FBmyo_IFN","FB_Postn_Thbs4",
  "FB_transition_Cd9", "FB_transition_Il1rapl1"
), zsc = T, cap = 2)


plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Ccl2",logsc=T,fr=F, um=T)




# cell cycle/proliferation dynamics

# Boxplot of G2M genes across time points
df <- colSums(sc1000_LAD_Cryo_Sham_Fibro@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_Fibro@counts)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
# plot normalized expressions across conditions
df1 <- data.frame(df[sample1], "Day01")
colnames(df1) <- c("G2M", "Days_post_MI")
df3 <- data.frame(df[sample3], "Day03")
colnames(df3) <- c("G2M", "Days_post_MI")
df7 <- data.frame(df[sample7], "Day07")
colnames(df7) <- c("G2M", "Days_post_MI")
df28 <- data.frame(df[sample28], "Day28")
colnames(df28) <- c("G2M", "Days_post_MI")
df56 <- data.frame(df[sample56], "Day56")
colnames(df56) <- c("G2M", "Days_post_MI")
dfSh <- data.frame(df[sampleSh], "Sham")
colnames(dfSh) <- c("G2M", "Days_post_MI")

df <- rbind(df1, df3, df7, df28, df56, dfSh)

df$Days_post_MI <- as.factor(df$Days_post_MI)

# box plot
# log2 scale
p<-ggplot(df, aes(x=Days_post_MI, y=log2(G2M), fill=Days_post_MI)) + geom_boxplot()+ 
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
# linear
p<-ggplot(df, aes(x=Days_post_MI, y=G2M, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 25)+
  stat_compare_means(method = "anova", label.y = 24)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))


# plot cell cycle scores
S_score   <- colMeans(sc1000_LAD_Cryo_Sham_Fibro@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_Fibro@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_Fibro@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_Fibro@ndata)),])

plotfeatmap(sc1000_LAD_Cryo_Sham_Fibro, G2M_score, "G2M score", logsc=T, um=T, cex=0.2)
plotfeatmap(sc1000_LAD_Cryo_Sham_Fibro, S_score, "S phase score", logsc=T, um=T, cex=0.2)







######




















# Transition probabilities
expData  <- getExpData(sc1000_LAD_Cryo_Sham_Fibro)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesFB, bmethod = "harmony")
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
# update cl after assigning new cluster identities
cl$partition <- sc1000_LAD_Cryo_Sham_Fibro@cpart

head(cl)
probs_Fb <-transitionProbs(res,cl,pvalue=0.01)
plotTrProbs(sc1000_LAD_Cryo_Sham_Fibro,probs_Fb,tp=.5,prthr=0,cthr=0,fr=F, um = T, cex=0.2)

# select clusters for pseudo-temporal plot (Slingshot)
require(FateID)
require(SingleCellExperiment)
require(slingshot)
require(DelayedMatrixStats)
set <- c(2,9,5) # from quiescent to D1 active
set <- c(2,9,5,3) # from quiescent to D1 active, then to D3 (not working)
set <- c(2,9,5,8) # from quiescent to D1 active, then to D3
set <- c(2,9,5,8,10) # from quiescent to D1 active, then to D3, then to D7
set <- c(2,9,5,8,10,9) # from quiescent to D1 active, then to D3, then to D7, then to transition, looping back to quiescent state
set <- c(2,9,5,8,15,10,9) # from quiescent to D1 active, then to D3, with cycling population, then to D7, then to transition, looping back to quiescent state (only works for umap)
set <- c(2,9,7,3,8,15) # from quiescent to D3 active
set <- c(2,9,7,3,8,15,10) # from quiescent to D3, then to D7 (not working)
set <- c(2,9,4,10) # from quiescent to D7 active
set <- c(10,4,9,2) # from D7 active back to quiescent
set <- c(5,3,15,10,4,9) # across active (myo)FBs
set <- c(9,7,5,3,15,10,4) # across active (myo)FBs (not working)

pt <- pseudoTime(sc1000_LAD_Cryo_Sham_Fibro,m="tsne",set=set)
# pt <- pseudoTime(sc1000_LAD_Cryo_Sham_Fibro,m="umap",set=set)
p1 <- plotPT(pt,sc1000_LAD_Cryo_Sham_Fibro,clusters=FALSE)
p2 <- plotPT(pt,sc1000_LAD_Cryo_Sham_Fibro,clusters=T)
p1 + p2

# plot gene module heat map
fs <- extractCounts(sc1000_LAD_Cryo_Sham_Fibro,minexpr=5,minnumber=5,pt=pt)
s1d   <- getsom(fs,nb=50,alpha=1)
ps    <- procsom(s1d,corthr=.85,minsom=0)
part  <- pt$part
ord   <- pt$ord
plotheatmap(ps$all.z, xpart=part[ord], xcol=sc1000_LAD_Cryo_Sham_Fibro@fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)
# gene list in each module:
module_genelist <- ps$nodes
names(module_genelist)[module_genelist%in%c(9)] # input module number
modules_dediff <- names(module_genelist)[module_genelist%in%c(9,12,14,15)]
# plotting only gene modules 1-9:
plotheatmap(ps$all.z[names(module_genelist)[module_genelist%in%c(9:13)],], xpart=part[ord], xcol=sc1000_LAD_Cryo_Sham_Fibro@fcol, ypart=ps$nodes[module_genelist%in%c(9:13)], xgrid=FALSE, ygrid=TRUE, xlab=TRUE)

# plot selected gene expression
plotexpression(fs,y=part,g="Erc2",n=ord,col=sc1000_LAD_Cryo_Sham_Fibro@fcol,cex=1,alpha=1)








# plotting cellular dynamics

# myofibroblast population dynamics (defined as clusters that are Postn+ and/or Acta2+; % of all fibroblasts)
myoFb_cells <- names(sc1000_LAD_Cryo_Sham_Fibro@cpart)[sc1000_LAD_Cryo_Sham_Fibro@cpart%in%c(2,6,7,8,9,15)]
sample1 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
v1 <- length(intersect(sample1, myoFb_cells)) / length(sample1) *100
v3 <- length(intersect(sample3, myoFb_cells)) / length(sample3) *100
v7 <- length(intersect(sample7, myoFb_cells)) / length(sample7) *100
v28 <- length(intersect(sample28, myoFb_cells)) / length(sample28) *100
v56 <- length(intersect(sample56, myoFb_cells)) / length(sample56) *100
vSh <- length(intersect(sampleSh, myoFb_cells)) / length(sampleSh) *100
a<-c(vSh,v1,v3,v7,v28,v56)
names(a)<-timepoints
barplot(a, xlab = "Days post-MI", ylab = "percentage of myofibroblasts",main = "% of myofibroblasts in fibroblast compartment", col = c("#CC33FF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000CC"))





# Dynamics of gene expression
# gene expression in a certain condition only
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
sample1 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
#sample <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_MI_D028|NM_Ad_LA_D56|NM_Ad_LA_D28|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Fibro@ndata))]
plotexpmap(sc1000_LAD_Cryo_Sham_Fibro,"Il10",logsc=T,fr=F, um=T, cex=1, cells = sample1)

# subset expression of gene of interest
df <- sc1000_LAD_Cryo_Sham_Fibro@ndata["Ccl19",]
# summed gene expression (normalized) in each time point
size1 <- length(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))])
size3 <- length(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))])
size7 <- length(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))])
size28 <- length(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))])
size56 <- length(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))])
sizeSh <- length(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))])

v1 <- sum(df[intersect(sample1, names(df))]) / size1
v3 <- sum(df[intersect(sample3, names(df))]) / size3
v7 <- sum(df[intersect(sample7, names(df))]) / size7
v28 <- sum(df[intersect(sample28, names(df))]) / size28
v56 <- sum(df[intersect(sample56, names(df))]) / size56
vSh <- sum(df[intersect(sampleSh, names(df))]) / sizeSh
a<-c(vSh,v1,v3,v7,v28,v56)
names(a)<-timepoints
barplot(a, xlab = "Days post-MI", ylab = "normalized gene expression",main = "Ccl19 expression in Fibroblasts", col = c("#CC33FF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000CC"))




