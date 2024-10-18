# immune cells in main cluster: 22,6,20,8,18,5,26,7,10,19,28,16,30,25,11,9
# subcluster Immune cells

Immune <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)[sc1000_LAD_Cryo_Sham_harmony@cpart%in%c(22,6,20,8,18,5,26,7,10,19,28,16,30,25,11,9,31)]
sc1000_LAD_Cryo_Sham_Immune1 <-SCseq(sc1000_LAD_Cryo_Sham_harmony@expdata[,Immune])
sc1000_LAD_Cryo_Sham_Immune1<-filterdata(sc1000_LAD_Cryo_Sham_Immune1,mintotal=1000, FGenes=rownames(sc1000_LAD_Cryo_Sham_Immune1@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(sc1000_LAD_Cryo_Sham_Immune1@expdata))])
expData  <- getExpData(sc1000_LAD_Cryo_Sham_Immune1)
S_score   <- colMeans(sc1000_LAD_Cryo_Sham_Immune1@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_Immune1@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_Immune1@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_Immune1@ndata)),])
regVar <- data.frame(S_score=S_score, G2M_score=G2M_score)
rownames(regVar) <- colnames(expData)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, regVar = regVar)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham_Immune1@expdata)[sc1000_LAD_Cryo_Sham_Immune1@genes,colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)]
noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=10)
sc1000_LAD_Cryo_Sham_Immune1 <- updateSC(sc1000_LAD_Cryo_Sham_Immune1,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham_Immune1 <- comptsne(sc1000_LAD_Cryo_Sham_Immune1)
sc1000_LAD_Cryo_Sham_Immune1 <- compumap(sc1000_LAD_Cryo_Sham_Immune1, n_neighbors = 50, min_dist = 0.8)
plotmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex = 0.2)

# plot cell cycle scores
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, "G2M score", logsc=T, um=T, cex=0.2)
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, S_score, "S phase score", logsc=T, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, cc_genes$s, "S phase score", um=T, logsc=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, cc_genes$g2m, "G2M phase score", um=T, logsc=T, cex=0.2)

sample1 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]

plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, logsc=T, um=T, cex=1, cells = sampleSh, n = "G2M score Sham")
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, logsc=T, um=T, cex=1, cells = sample1, n = "G2M score Day 1")
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, logsc=T, um=T, cex=1, cells = sample3, n = "G2M score Day 3")
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, logsc=T, um=T, cex=1, cells = sample7, n = "G2M score Day 7")
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, logsc=T, um=T, cex=1, cells = sample28, n = "G2M score Day 28")
plotfeatmap(sc1000_LAD_Cryo_Sham_Immune1, G2M_score, logsc=T, um=T, cex=1, cells = sample56, n = "G2M score Day 56")

# Violin plot G2M genes across clusters
violinMarkerPlot(cc_genes$g2m,sc1000_LAD_Cryo_Sham_Immune1,set=c(8,7,6,13,2,1,20,22), ti = "G2M genes")

# Plot violin plot of G2M genes across time points
Macrophages <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart%in%c(8,7,6,13,2,1,20,22)]
df <- colSums(sc1000_LAD_Cryo_Sham_Immune1@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_Immune1@counts)
sample1 <- intersect(sample1, Macrophages)
sample3 <- intersect(sample3, Macrophages)
sample7 <- intersect(sample7, Macrophages)
sample28 <- intersect(sample28, Macrophages)
sample56 <- intersect(sample56, Macrophages)
sampleSh <- intersect(sampleSh, Macrophages)
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


# violin plot
require(ggpubr)
p<-ggplot(df, aes(x=Days_post_MI, y=log2(G2M), fill=Days_post_MI)) + geom_violin()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
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




# Boxplot of S phase genes across time points
Macrophages <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart%in%c(8,7,6,13,2,1,20,22)]
df <- colSums(sc1000_LAD_Cryo_Sham_Immune1@ndata[cc_genes$s,]) * median(sc1000_LAD_Cryo_Sham_Immune1@counts)
sample1 <- intersect(sample1, Macrophages)
sample3 <- intersect(sample3, Macrophages)
sample7 <- intersect(sample7, Macrophages)
sample28 <- intersect(sample28, Macrophages)
sample56 <- intersect(sample56, Macrophages)
sampleSh <- intersect(sampleSh, Macrophages)
# plot normalized expressions across conditions
df1 <- data.frame(df[sample1], "Day01")
colnames(df1) <- c("S", "Days_post_MI")
df3 <- data.frame(df[sample3], "Day03")
colnames(df3) <- c("S", "Days_post_MI")
df7 <- data.frame(df[sample7], "Day07")
colnames(df7) <- c("S", "Days_post_MI")
df28 <- data.frame(df[sample28], "Day28")
colnames(df28) <- c("S", "Days_post_MI")
df56 <- data.frame(df[sample56], "Day56")
colnames(df56) <- c("S", "Days_post_MI")
dfSh <- data.frame(df[sampleSh], "Sham")
colnames(dfSh) <- c("S", "Days_post_MI")

df <- rbind(df1, df3, df7, df28, df56, dfSh)

df$Days_post_MI <- as.factor(df$Days_post_MI)

# box plot
# log2 scale
p<-ggplot(df, aes(x=Days_post_MI, y=log2(S), fill=Days_post_MI)) + geom_boxplot()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
# linear
p<-ggplot(df, aes(x=Days_post_MI, y=S, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 7.5)+
  stat_compare_means(method = "anova", label.y = 6)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))






# Boxplot of G2M genes across time points, only among resident MPs
Macrophages_resident <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart%in%c(2,1,20,22)]
df <- colSums(sc1000_LAD_Cryo_Sham_Immune1@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_Immune1@counts)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample1 <- intersect(sample1, Macrophages_resident)
sample3 <- intersect(sample3, Macrophages_resident)
sample7 <- intersect(sample7, Macrophages_resident)
sample28 <- intersect(sample28, Macrophages_resident)
sample56 <- intersect(sample56, Macrophages_resident)
sampleSh <- intersect(sampleSh, Macrophages_resident)
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



# Boxplot of G2M genes across time points, only among ciculating MPs
Macrophages_circulating <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart%in%c(6,7,8,13)]
df <- colSums(sc1000_LAD_Cryo_Sham_Immune1@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_Immune1@counts)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata))]
sample1 <- intersect(sample1, Macrophages_circulating)
sample3 <- intersect(sample3, Macrophages_circulating)
sample7 <- intersect(sample7, Macrophages_circulating)
sample28 <- intersect(sample28, Macrophages_circulating)
sample56 <- intersect(sample56, Macrophages_circulating)
sampleSh <- intersect(sampleSh, Macrophages_circulating)
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













# Transition probabilities, Slingshot pseudotime and gene modules
probs_Immune <-transitionProbs(res,cl,pvalue=0.01)
plotTrProbs(sc1000_LAD_Cryo_Sham_Immune1,probs_Immune,tp=.5,prthr=0,cthr=0,fr=F, um = T, cex=0.05)

# select clusters for pseudo-temporal plot (Slingshot)
require(FateID)
require(SingleCellExperiment)
require(slingshot)
require(DelayedMatrixStats)
set <- c(8,7,6,13,2,1) # Whole Macrophage lineage
set <- c(7,6,13,2,1) # Macrophage lineage
set <- c(8,7,6,13) # Early Macrophage lineage
pt <- pseudoTime(sc1000_LAD_Cryo_Sham_Immune1,m="umap",set=set)
plotPT(pt,sc1000_LAD_Cryo_Sham_Immune1,clusters=FALSE)
plotPT(pt,sc1000_LAD_Cryo_Sham_Immune1,clusters=T)

# plot gene module heat map
fs <- extractCounts(sc1000_LAD_Cryo_Sham_Immune1,minexpr=5,minnumber=20,pt=pt)
s1d   <- getsom(fs, nb=50, alpha=1)
ps    <- procsom(s1d,corthr=.85,minsom=3)
part  <- pt$part
ord   <- pt$ord
plotheatmap(ps$all.z, xpart=part[ord], xcol=sc1000_LAD_Cryo_Sham_Immune1@fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)
# plot Z-scores
plotheatmap(ps$all.b,xpart=part[ord],xcol=sc1000_LAD_Cryo_Sham_Immune1@fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
# gene list in each module:
module_genelist <- ps$nodes
names(module_genelist)[module_genelist%in%c(9)] # input module number
modules_macrotransition <- names(module_genelist)[module_genelist%in%c(9)]
# plotting only gene module 9:
plotheatmap(ps$all.z[names(module_genelist)[module_genelist%in%c(9)],], xpart=part[ord], xcol=sc1000_LAD_Cryo_Sham_Immune1@fcol, ypart=ps$nodes[module_genelist%in%c(9)], xgrid=FALSE, ygrid=TRUE, xlab=TRUE)
plotexpression(fs,part[ord],g,n=ord,col=sc1000_LAD_Cryo_Sham_Immune1@fcol,name="Node 9",cluster=FALSE,alpha=.5,types=NULL)

plotexpressionProfile(fs,y=part,g=modules_macrotransition,n=ord,alpha=1,col=rainbow(length(modules_macrotransition)),lwd = 2, cluster = T,logsc = T)

# find certain gene in module genelist
module_genelist["Gpnmb"]

# plot selected gene expression
plotexpression(fs,y=part,g="Gpnmb",n=ord,col=sc1000_LAD_Cryo_Sham_Immune1@fcol,cex=1,alpha=1)
plotexpression(fs,y=part,g="Cd274",n=ord,col=sc1000_LAD_Cryo_Sham_Immune1@fcol,cex=1,alpha=1)
plotexpression(fs,y=part,g="Mmp9",n=ord,col=sc1000_LAD_Cryo_Sham_Immune1@fcol,cex=1,alpha=1)









# change cluster colors
library(ggplot2)
library(paletteer)
sc1000_LAD_Cryo_Sham_Immune1_fcol_original <- sc1000_LAD_Cryo_Sham_Immune1@fcol
# Neutro: 3, 4
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(3,4)] <- c("#FF7043", "#BF360C")
# Mast cells: 24
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(24)] <- "#4E342E"
# Erythrocytes: 23
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(23)] <- "#FFCCBC"
# Mono-Macro: 8,7,6,13,2,1,20,22
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(8,7,6,13,2,1,20,22)] <- as.character(paletteer_c("grDevices::PuOr", 30))[c(4,8,10,12,23,25,27,30)]
# DCs: 17,5,21,19,11
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(17,5,21,11,19)] <- c(as.character(paletteer_d("ggsci::teal_material", 10))[c(3,6,7,9)], "#004D40")
# B cells: 10
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(10)] <- "#D81B60"
# T cells: 18,9,12
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(18,9,12)] <- as.character(paletteer_d("ggsci::purple_material", 10))[c(5,7,9)]
# NK cells: 14; # ILC2: 15
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(14,15)] <- c("#3949AB", "#5C6BC0")
# Fibrocytes: 16
sc1000_LAD_Cryo_Sham_Immune1@fcol[c(16)] <- "#B0BEC5"

plotmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex=0.2)

plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Fn1",logsc=T,fr=F, um=T, cex=0.2)
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,c(13),pvalue=0.01)


# plot umap with cluster names annotated
require(paletteer)
cluster_ann <- sc1000_LAD_Cryo_Sham_Immune1@cpart
cluster_ann[cluster_ann%in%1] <- "MP_Timd4_Lyve1"
cluster_ann[cluster_ann%in%2] <- "MP_MHCII_Cx3cr1"
cluster_ann[cluster_ann%in%3] <- "Neutrophil_1" 
cluster_ann[cluster_ann%in%4] <- "Neutrophil_2" 
cluster_ann[cluster_ann%in%5] <- "DCc2b_moDC" 
cluster_ann[cluster_ann%in%6] <- "MP_Spp1"
cluster_ann[cluster_ann%in%7] <- "MP_Ly6C_mid" 
cluster_ann[cluster_ann%in%8] <- "MP_Ly6C_high"
cluster_ann[cluster_ann%in%9] <- "T_naive"
cluster_ann[cluster_ann%in%10] <- "B_cells"
cluster_ann[cluster_ann%in%11] <- "DCc1" 
cluster_ann[cluster_ann%in%12] <- "T_Eff"
cluster_ann[cluster_ann%in%13] <- "MP_Trem2" 
cluster_ann[cluster_ann%in%14] <- "NK"
cluster_ann[cluster_ann%in%15] <- "ILC2"
cluster_ann[cluster_ann%in%16] <- "Fibrocytes" 
cluster_ann[cluster_ann%in%17] <- "DC_Ifitm1" 
cluster_ann[cluster_ann%in%18] <- "T_gd" 
cluster_ann[cluster_ann%in%19] <- "DC_Ccr7" 
cluster_ann[cluster_ann%in%20] <- "MP_Retnla"
cluster_ann[cluster_ann%in%21] <- "DCp_Cd8" 
cluster_ann[cluster_ann%in%22] <- "MP_Cxcl13" 
cluster_ann[cluster_ann%in%23] <- "Erythrocyte" 
cluster_ann[cluster_ann%in%24] <- "Mast_cells"
cluster_ann_immune <- cluster_ann
immune_manual_col <- c("#D81B60", c(as.character(paletteer_d("ggsci::teal_material", 10))[c(6,7,3,9)], "#004D40"), "#FFCCBC", "#B0BEC5", "#5C6BC0", as.character(paletteer_c("grDevices::PuOr", 30))[c(30,23,27,25,12)], "#4E342E", as.character(paletteer_c("grDevices::PuOr", 30))[c(8,4,10)], c("#FF7043", "#BF360C"), "#3949AB", as.character(paletteer_d("ggsci::purple_material", 10))[c(9,5,7)] 
)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann_immune, um=T, cex=0.2, samples_col = 
                 immune_manual_col
)



plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann_immune, um=T, cex=0.2, samples_col = sc1000_LAD_Cryo_Sham_Immune1@fcol,
               , leg=T)


immune_populations <- cluster_ann_immune

# matching cell type names with colors
sc1000_LAD_Cryo_Sham_Immune1@fcol
head(sc1000_LAD_Cryo_Sham_Immune1@cpart)

Immune_fcol <- sc1000_LAD_Cryo_Sham_Immune1@fcol[sc1000_LAD_Cryo_Sham_Immune1@cpart]
names(Immune_fcol) <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)
head(Immune_fcol)

head(immune_populations)
identical(names(immune_populations), names(Immune_fcol))

Immune_clusters_fcol <- data.frame(clusters=immune_populations, colors=Immune_fcol)

saveRDS(Immune_clusters_fcol, "df_cluster_colors_Immune.rds")





# plot macrophages only
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann_immune, um=T, cex=0.2, samples_col = 
                 c("#D81B60", 
                   c(as.character(paletteer_d("ggsci::teal_material", 10))[c(6,7,3,9)], "#004D40"), 
                   "#FFCCBC", "#B0BEC5", "#5C6BC0", as.character(paletteer_c("grDevices::PuOr", 30))[c(30,23,27,25,12)], "#4E342E", as.character(paletteer_c("grDevices::PuOr", 30))[c(8,4,10)], c("#FF7043", "#BF360C"), "#3949AB", as.character(paletteer_d("ggsci::purple_material", 10))[c(9,5,7)] 
                 ),
               subset = c("MP_Cxcl13" ,"MP_MHCII_Cx3cr1", "MP_Retnla","MP_Timd4_Lyve1", "MP_Trem2", "MP_Ly6C_mid", "MP_Ly6C_high", "MP_Spp1")
               , leg=T)

plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann_immune, um=T, cex=0.2, samples_col = sc1000_LAD_Cryo_Sham_Immune1@fcol,
               subset = c("22 MP_Cxcl13" ,"02 MP_MHCII_Cx3cr1", "20 MP_Retnla","01 MP_Timd4_Lyve1", "13 MP_Trem2", "07 MP_Ly6C_mid", "08 MP_Ly6C_high", "06 MP_Spp1")
               , leg=F)


# plot DCs only
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann, um=T, cex=0.2, samples_col = 
                 c("#D81B60", c(as.character(paletteer_d("ggsci::teal_material", 10))[c(6,7,3,9)], "#004D40"), "#FFCCBC", "#B0BEC5", "#5C6BC0", as.character(paletteer_c("grDevices::PuOr", 30))[c(30,23,27,25,12)], "#4E342E", as.character(paletteer_c("grDevices::PuOr", 30))[c(8,4,10)], c("#FF7043", "#BF360C"), "#3949AB", as.character(paletteer_d("ggsci::purple_material", 10))[c(9,5,7)] 
                 ),
               subset = c("cDC1","cDC2a", "cDC2b", "DC_Ccr7", "DCp_Cd8" ))

# plot Lymphocytes only
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann, um=T, cex=0.2, samples_col = 
                 c("#D81B60", c(as.character(paletteer_d("ggsci::teal_material", 10))[c(6,7,3,9)], "#004D40"), "#FFCCBC", "#B0BEC5", "#5C6BC0", as.character(paletteer_c("grDevices::PuOr", 30))[c(30,23,27,25,12)], "#4E342E", as.character(paletteer_c("grDevices::PuOr", 30))[c(8,4,10)], c("#FF7043", "#BF360C"), "#3949AB", as.character(paletteer_d("ggsci::purple_material", 10))[c(9,5,7)] 
                 ),
               subset = c("B_cells", "ILC2", "NK", "T_Eff", "T_gd", "T_naive"))

# plot Neutrophils only
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, cluster_ann, um=T, cex=0.2, samples_col = 
                 c("#D81B60", c(as.character(paletteer_d("ggsci::teal_material", 10))[c(6,7,3,9)], "#004D40"), "#FFCCBC", "#B0BEC5", "#5C6BC0", as.character(paletteer_c("grDevices::PuOr", 30))[c(30,23,27,25,12)], "#4E342E", as.character(paletteer_c("grDevices::PuOr", 30))[c(8,4,10)], c("#FF7043", "#BF360C"), "#3949AB", as.character(paletteer_d("ggsci::purple_material", 10))[c(9,5,7)] 
                 ),
               subset = c("Neutrophil_1", "Neutrophil_2"))


# gene expression heatmap

plotmarkergenes(sc1000_LAD_Cryo_Sham_Immune1,c("Pcna",cc_genes$s[c(1:10)], "Mki67",cc_genes$g2mc[c(1:10)]),cl=c(8,7,6,13,2,1,20,22),samples=types,order.cells=F, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"), fontsize = 10)
fractDotPlot(sc1000_LAD_Cryo_Sham_Immune1, c(cc_genes$s, cc_genes$g2m), cl=c(8,7,6,13,2,1,20,22), zsc=TRUE)

dg8 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,8,pvalue=.01)
dg8 <- dg8$dg
genes8 <- dg8[dg8$fc>1, c("fc", "padj")]
genes8 <- head(rownames(genes8)[order(genes8$fc, decreasing = T)],5)


dg7 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,7,pvalue=.01)
dg7 <- dg7$dg
genes7 <- dg7[dg7$fc>1, c("fc", "padj")]
genes7 <- head(rownames(genes7)[order(genes7$fc, decreasing = T)],5)

dg6 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,6,pvalue=.01)
dg6 <- dg6$dg
genes6 <- dg6[dg6$fc>1, c("fc", "padj")]
genes6 <- head(rownames(genes6)[order(genes6$fc, decreasing = T)],5)

dg13 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,13,pvalue=.01)
dg13 <- dg13$dg
genes13 <- dg13[dg13$fc>1, c("fc", "padj")]
genes13 <- head(rownames(genes13)[order(genes13$fc, decreasing = T)],5)

dg2 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,2,pvalue=.01)
dg2 <- dg2$dg
genes2 <- dg2[dg2$fc>1, c("fc", "padj")]
genes2 <- head(rownames(genes2)[order(genes2$fc, decreasing = T)],5)

dg1 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,1,pvalue=.01)
dg1 <- dg1$dg
genes1 <- dg1[dg1$fc>1, c("fc", "padj")]
genes1 <- head(rownames(genes1)[order(genes1$fc, decreasing = T)],5)

dg20 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,20,pvalue=.01)
dg20 <- dg20$dg
genes20 <- dg20[dg20$fc>1, c("fc", "padj")]
genes20 <- head(rownames(genes20)[order(genes20$fc, decreasing = T)],5)

dg22 <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,22,pvalue=.01)
dg22 <- dg22$dg
genes22 <- dg22[dg22$fc>1, c("fc", "padj")]
genes22 <- head(rownames(genes22)[order(genes22$fc, decreasing = T)],5)

genes <- unique(c(genes8, genes7, genes6, genes13, genes2, genes1, genes20, genes22))
plotmarkergenes(sc1000_LAD_Cryo_Sham_Immune1,genes,cl=c(8,7,6,13,2,1,20,22),samples=types,order.cells=F, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"), fontsize = 4)


# cytokine profile
# list1: Tnf, Il6, Il1b, Mmp9, Mmp13, Il10, Tgfb1, Il10ra, Ccl7, Ccl8
# list2: Il1a, Il1b, Il6, Il8, Il12, Il17, Il18, Il20, Il23, Il33, Tnf, Ifng, Tgfb1, Lif, Csf2, Osm, Il4, Il10, Il11, Il13, Lta
cytokine_list1 <- c("Tnf", "Il6", "Il1b", "Mmp9", "Mmp13", "Il10", "Tgfb1", "Il10ra", "Ccl7", "Ccl8")
cytokine_list2 <- unique(c("Il1a", "Il1b", "Il6", "Il12a", "Il17a", "Il18", "Il20", "Il23a", "Il33", "Tnf", "Ifng", "Tgfb1", "Lif", "Csf2", "Osm", "Il4", "Il10", "Il11", "Il13", "Lta"))
cytokines_proinf <- unique(c("Il1a", "Il1b", "Il6", "Il12a", "Il17a", "Il18", "Il20", "Il23a", "Il33", "Tnf", "Ifng", "Tgfb1", "Lif", "Csf2", "Osm"))
cytokines_antiinf<- c("Il4", "Il10", "Il11", "Il13", "Lta")

plotmarkergenes(sc1000_LAD_Cryo_Sham_Immune1,cytokine_list1,cl=c(8,7,6,13,2,1,20,22),samples=types,order.cells=F, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"), fontsize = 6)
plotmarkergenes(sc1000_LAD_Cryo_Sham_Immune1,cytokine_list2,cl=c(8,7,6,13,2,1,20,22),samples=types,order.cells=F, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"), fontsize = 6)

plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, cytokines_proinf, n = "Proinflammatory cytokine score", um=T, cex=0.2, logsc=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, cytokines_antiinf, n = "Anti-inflammatory cytokine score", um=T, cex=1, logsc=T)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, "Il5ra", um=T, cex=0.2, logsc=T)

plotexpmap(sc1000_LAD_Cryo_Sham_Immune1, cytokines_proinf, n = "Proinflammatory cytokine score", um=T, cex=1, logsc=T, cells=sample7)





# visualize different conditions
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
typesC <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata), 7,8)
head(typesC)
unique(typesC)
typesC[typesC%in%"LA"] <- "LAD"
typesC[typesC%in%"MI"] <- "Cryoablation"
typesC[typesC%in%"Sh"] <- "Sham"
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesC, um=T, cex=0.05, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesC, um=T, subset = "LAD")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesC, um=T, subset = "Cryoablation")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesC, um=T, subset = "Sham", samples_col  = rep("#0A1722", length(typesC[typesC%in%"Sham"])))

# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
typesT <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Immune1@ndata), 7,12)

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
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, cex = 0.2, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "D01")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "D03")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "D07")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "D28")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "Sham")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "LA_D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "MI_D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = "Sh_D56")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = c("Sh_D01","Sh_D03","Sh_D07","Sh_D28", "Sh_D56"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Immune1, typesT, um=T, subset = c("D01", "D03", "Sham"), samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
# immune marker genes
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"S100a8",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Itgam",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Itgax",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Lyve1",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cx3cr1",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Retnla",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Itgae",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cd209a",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cd3e",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Il7r",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Il2ra",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Gata3",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Gzma",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cd19",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Mrc1",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Itgal",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Spp1",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Ccr2",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cd4",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cd8a",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Foxp3",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Trdc",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Hbb-bt",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cpa3",logsc=T,fr=F, um=T, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Trem2",logsc=T,fr=F, um=T, cex=0.2)


A <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart %in% c(2)]
B <- names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart %in% c(8)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_Immune1,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="D56_Shared_EC",Bname="D56_LAD_EC",show_names=TRUE,padj=TRUE)

plotexpmap(sc1000_LAD_Cryo_Sham_Immune1,"Cdh5",logsc=T,fr=F, um=T)

dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,c(21),pvalue=0.01)







# dotplots of selected genes
# pool clusters into defined sub-populations
plotmap(sc1000_LAD_Cryo_Sham_Immune1, um=T, cex=0.5)
Immune_populations <- sc1000_LAD_Cryo_Sham_Immune1@cpart
Immune_populations[Immune_populations %in% c(3,4)] <- "Neutrophils"
Immune_populations[Immune_populations %in% c(8,7)] <- "Ly6C_high/mid MP"
Immune_populations[Immune_populations %in% c(6)] <- "Spp1+ MP"
Immune_populations[Immune_populations %in% c(13)] <- "Trem2+ MP"
Immune_populations[Immune_populations %in% c(2)] <- "Cx3cr1_high MP"
Immune_populations[Immune_populations %in% c(1)] <- "Timd4+ MP"
Immune_populations[Immune_populations %in% c(20)] <- "Retnla+ MP"
Immune_populations[Immune_populations %in% c(22)] <- "Cxcl13+ MP"
Immune_populations[Immune_populations %in% c(11)] <- "cDC1"
Immune_populations[Immune_populations %in% c(17,5)] <- "cDC2/moDC"
Immune_populations[Immune_populations %in% c(19)] <- "CCR7 DC"
Immune_populations[Immune_populations %in% c(21)] <- "pDC"
Immune_populations[Immune_populations %in% c(24)] <- "MC"
Immune_populations[Immune_populations %in% c(15)] <- "ILC2"
Immune_populations[Immune_populations %in% c(9,12)] <- "T cells"
Immune_populations[Immune_populations %in% c(18)] <- "gdT cells"
Immune_populations[Immune_populations %in% c(14)] <- "NK cells"
Immune_populations[Immune_populations %in% c(10)] <- "B cells"
Immune_populations[Immune_populations %in% c(16)] <- "Fibrocytes"
Immune_populations[Immune_populations %in% c(23)] <- "Erythrocytes"


table(Immune_populations)
# interested genes
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Immune1,c(2,3),pvalue=0.01)

genes <- c(
  "S100a8", "Ly6c2", "Ccr2", "Itgam", "Spp1", "Trem2", "Mrc1",
  "H2-Ab1", "Cx3cr1", "Timd4", "Lyve1", "Retnla", "Cxcl13",
  "Itgax", "Itgae", "Cd209a", "Mgl2", "Ccr7", "Siglech",
  "Cpa3",
  "Il7r", "Il2ra", "Gata3", "Cd3e", "Cd4", "Cd8a", "Foxp3", "Trdc", 
  "Gzma",
  "Cd19",
  "Hbb-bt"
)

fractDotPlot(sc1000_LAD_Cryo_Sham_Immune1, genes, samples = Immune_populations, logscale=TRUE)
fractDotPlot(sc1000_LAD_Cryo_Sham_Immune1, genes, samples = Immune_populations, subset=c("Neutrophils", "Ly6C_high/mid MP", "Spp1+ MP", "Trem2+ MP",
                                                                                         "Cx3cr1_high MP", "Timd4+ MP", "Retnla+ MP", "Cxcl13+ MP",
                                                                                         "cDC1", "cDC2/moDC", "CCR7 DC", "pDC", "MC", "ILC2", "T cells",
                                                                                         "gdT cells", "NK cells", "B cells", "Fibrocytes", "Erythrocytes"), zsc = T, cap = 2)

# simplify genes and cell type lists
genes1 <- c(
  "S100a8", "Ly6c2", "Ccr2", "Itgam",  
  "H2-Ab1", "Cx3cr1", "Timd4", "Lyve1",
  "Itgae", "Cd209a", "Mgl2", "Ccr7", "Siglech",
  "Gata3", "Cd3e", "Cd4", "Cd8a", "Foxp3",  "Gzma",  "Cd19"
)
fractDotPlot(sc1000_LAD_Cryo_Sham_Immune1, genes1, samples = Immune_populations, subset=c("Neutrophils", "Ly6C_high/mid MP",
                                                                                          "Timd4+ MP",
                                                                                          "cDC1", "cDC2/moDC", "CCR7 DC", "pDC", "ILC2", "T cells")













