# subcluster SwC populations
Neural <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)[sc1000_LAD_Cryo_Sham_harmony@cpart%in%c(29)]
sc1000_LAD_Cryo_Sham_neural <-SCseq(sc1000_LAD_Cryo_Sham_harmony@expdata[,Neural])
sc1000_LAD_Cryo_Sham_neural<-filterdata(sc1000_LAD_Cryo_Sham_neural,mintotal=1000, FGenes=rownames(sc1000_LAD_Cryo_Sham_neural@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(sc1000_LAD_Cryo_Sham_neural@expdata))])
expData  <- getExpData(sc1000_LAD_Cryo_Sham_neural)
# with batch removal (harmony)
batches <- Neural
names(batches) <- Neural
batch1 <- batches[grep("D01|D03|D07|D28|MI_D56",batches)]
batch2 <- batches[grep("LA_D56|Sh_D56",batches)]
batch1 <- replace(batch1,,"b1")
batch2 <- replace(batch2,,"b2")
batchesNeural <- c(batch1, batch2)

S_score   <- colMeans(sc1000_LAD_Cryo_Sham_neural@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_neural@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_neural@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_neural@ndata)),])
regVar <- data.frame(S_score=S_score, G2M_score=G2M_score)
rownames(regVar) <- colnames(expData)

res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesNeural, bmethod = "harmony", regVar = regVar)
# no batch effect removal
# res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham_neural@expdata)[sc1000_LAD_Cryo_Sham_neural@genes,colnames(sc1000_LAD_Cryo_Sham_neural@ndata)]
# noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=10)
sc1000_LAD_Cryo_Sham_neural <- updateSC(sc1000_LAD_Cryo_Sham_neural,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham_neural <- comptsne(sc1000_LAD_Cryo_Sham_neural)
sc1000_LAD_Cryo_Sham_neural <- compumap(sc1000_LAD_Cryo_Sham_neural, n_neighbors = 10, min_dist = 0.8)
#sc1000_LAD_Cryo_Sham_neural <- compumap(sc1000_LAD_Cryo_Sham_neural, n_neighbors = 10, min_dist = 0.5)
# sc1000_LAD_Cryo_Sham_neural <- compumap(sc1000_LAD_Cryo_Sham_neural)
plotmap(sc1000_LAD_Cryo_Sham_neural, um=T, cex = 1)
plotmap(sc1000_LAD_Cryo_Sham_neural, um=F, cex = 0.2)
table(sc1000_LAD_Cryo_Sham_neural@cpart)
# change cluster colors
library(ggplot2)
library(paletteer)
sc1000_LAD_Cryo_Sham_neural@fcol[c(3)] <- "#48C2BE"
sc1000_LAD_Cryo_Sham_neural@fcol[c(2)] <- "#5B7CBC"
sc1000_LAD_Cryo_Sham_neural@fcol[c(1)] <- "#804AA4"
sc1000_LAD_Cryo_Sham_neural@fcol[c(4)] <- "#4C8D56"
plotmap(sc1000_LAD_Cryo_Sham_neural, um=T, cex=1)

dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_neural,c(2),pvalue=0.5)



plotexpmap(sc1000_LAD_Cryo_Sham_neural,"Plp1",logsc=T,fr=F, um=T, cex=2)
plotexpmap(sc1000_LAD_Cryo_Sham_neural,"Prnp",logsc=T,fr=F, um=T, cex=2)
plotexpmap(sc1000_LAD_Cryo_Sham_neural,"Pdgfrb",logsc=T,fr=F, um=T, cex=2)
plotexpmap(sc1000_LAD_Cryo_Sham_neural,"Pdgfra",logsc=T,fr=F, um=T, cex=2)

# plot only a certain condition/timepoint
sampleD28 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("MI_D28|LA_D28",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sampleD7 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("MI_D07|LA_D07",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("Sh",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
plotexpmap(sc1000_LAD_Cryo_Sham_neural,"Wnt6",logsc=T,fr=F, um=T, cex=4, cells = sampleSh)


# plot umap with cluster names annotated

cluster_ann <- sc1000_LAD_Cryo_Sham_neural@cpart
cluster_ann[cluster_ann%in%1] <- "Schwann_quiescent"
cluster_ann[cluster_ann%in%2] <- "Schwann_Galectin"
cluster_ann[cluster_ann%in%3] <- "Schwann_IFN" 
cluster_ann[cluster_ann%in%4] <- "Schwann_metabolic" 

cluster_ann_neural <- cluster_ann

plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, cluster_ann_neural, um=T, cex=1, samples_col = 
                 sc1000_LAD_Cryo_Sham_neural@fcol
)



neural_populations <- cluster_ann_neural




# visualize different conditions
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_neural@ndata), 7,8)
head(types)
unique(types)
types[types%in%"LA"] <- "LAD"
types[types%in%"MI"] <- "Cryoablation"
types[types%in%"Sh"] <- "Sham"
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, cex=1, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "LAD")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "Cryoablation")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "Sham", samples_col  = rep("#0A1722", length(types[types%in%"Sham"])), cex=2)

# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_neural@ndata), 7,12)

plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = c("Sh_D01","Sh_D03","Sh_D07","Sh_D28", "Sh_D56"))

# merge all Sham time points into 1 type
types<- gsub("Sh_D01", "Sham", types)
types<- gsub("Sh_D03", "Sham", types)
# types<- gsub("Sh_D07", "Sham", types) contaminated with LAD cells
types<- gsub("Sh_D28", "Sham", types)
types<- gsub("Sh_D56", "Sham", types)

types<- gsub("LA_D01", "D01", types)
types<- gsub("LA_D03", "D03", types)
types<- gsub("LA_D07", "D07", types)
types<- gsub("LA_D28", "D28", types)
types<- gsub("LA_D56", "D56", types)

types<- gsub("MI_D01", "D01", types)
types<- gsub("MI_D03", "D03", types)
types<- gsub("MI_D07", "D07", types)
types<- gsub("MI_D28", "D28", types)
types<- gsub("MI_D56", "D56", types)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, cex = 1, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "D01")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "D03")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "D07")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "D28")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "Sham")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "LA_D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "MI_D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = "Sh_D56")

plotsymbolsmap(sc1000_LAD_Cryo_Sham_neural, types, um=T, subset = c("Sh_D01","Sh_D03","Sh_D07","Sh_D28", "Sh_D56"))



# dotplots of selected genes
# pool clusters into defined sub-populations
plotmap(sc1000_LAD_Cryo_Sham_neural, um=T, cex=0.5)
SwC_populations <- cluster_ann_neural
table(SwC_populations)

# interested genes
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_neural,c(2),pvalue=0.01)
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_neural,c(2),pvalue=1)
genes <- c("Plp1", "Prnp", "Pdgfrb", "Apoe", "Cd63", "Lgals3", "Mt2","Isg15", "Cxcl10", "Ifi27l2a",
           "Lars2", "Camk1d", "Cmss1", "Hexb",
           "Mki67", "Pcna"
)

fractDotPlot(sc1000_LAD_Cryo_Sham_neural, genes, samples = SwC_populations, subset=c(
  "Schwann_quiescent", "Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic"
), zsc = T, cap = 2)


plotexpmap(sc1000_LAD_Cryo_Sham_neural,"Dct",logsc=T,fr=F, um=T, cex=2)




# Boxplot of G2M genes across time points
df <- colSums(sc1000_LAD_Cryo_Sham_neural@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_neural@counts)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_neural@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_neural@ndata))]

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
p<-ggplot(df, aes(x=Days_post_MI, y=G2M, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 9)+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))

plotfeatmap(sc1000_LAD_Cryo_Sham_neural, G2M_score, logsc=T, um=T, cex=2, n = "G2M score")




















# Transition probabilities and Slingshot pseudotime
expData  <- getExpData(sc1000_LAD_Cryo_Sham_neural)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesFB, bmethod = "harmony")
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
# update cl after assigning new cluster identities
cl$partition <- sc1000_LAD_Cryo_Sham_neural@cpart

head(cl)
probs_DC <-transitionProbs(res,cl,pvalue=0.01)
plotTrProbs(sc1000_LAD_Cryo_Sham_neural,probs_DC,tp=.5,prthr=0,cthr=0,fr=F, um = T, cex=2)

# select clusters for pseudo-temporal plot (Slingshot)
require(FateID)
require(SingleCellExperiment)
require(slingshot)
require(DelayedMatrixStats)
set <- c(1,2,3) # 


pt <- pseudoTime(sc1000_LAD_Cryo_Sham_neural,m="umap",set=set)
plotPT(pt,sc1000_LAD_Cryo_Sham_neural,clusters=FALSE)
plotPT(pt,sc1000_LAD_Cryo_Sham_neural,clusters=T)

