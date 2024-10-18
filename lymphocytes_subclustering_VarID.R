# assign new T cells, NK and ILC clusters according to higher definition subclustering
nonB_Lymphocytes <- names(sc1000_LAD_Cryo_Sham_Immune@cpart)[sc1000_LAD_Cryo_Sham_Immune@cpart%in%c(8,13,15,16,20,24)]
#identical(nonB_Lymphocytes, names(sc1000_LAD_Cryo_Sham_Immune1@cpart)[sc1000_LAD_Cryo_Sham_Immune1@cpart%in%c(15,18,14,12,9)])
sc1000_LAD_Cryo_Sham_nonBlymphocytes1 <-SCseq(sc1000_LAD_Cryo_Sham_Immune@expdata[,nonB_Lymphocytes])
sc1000_LAD_Cryo_Sham_nonBlymphocytes1<-filterdata(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,mintotal=1000, FGenes=rownames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@expdata))])
expData  <- getExpData(sc1000_LAD_Cryo_Sham_nonBlymphocytes1)

# batch effect removal and cell cycle genes regression
S_score   <- colMeans(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)),])
regVar <- data.frame(S_score=S_score, G2M_score=G2M_score)
rownames(regVar) <- colnames(expData)

batches <- nonB_Lymphocytes
names(batches) <- nonB_Lymphocytes
batch1 <- batches[grep("D01|D03|D07|D28|MI_D56",batches)]
batch2 <- batches[grep("LA_D56|Sh_D56",batches)]
batch1 <- replace(batch1,,"b1")
batch2 <- replace(batch2,,"b2")
batchesLymph <- c(batch1, batch2)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesLymph, bmethod = "harmony", regVar=regVar)
# no batch effect removal
# res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@expdata)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@genes,colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)]
# noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=10)
sc1000_LAD_Cryo_Sham_nonBlymphocytes1 <- updateSC(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham_nonBlymphocytes1 <- comptsne(sc1000_LAD_Cryo_Sham_nonBlymphocytes1)
sc1000_LAD_Cryo_Sham_nonBlymphocytes1 <- compumap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, spread = 3, min_dist = 1)
# sc1000_LAD_Cryo_Sham_nonBlymphocytes <- compumap(sc1000_LAD_Cryo_Sham_nonBlymphocytes)
plotmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, um=T, cex = 0.5)
plotmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, um=F, cex = 0.2)
table(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)



# change cluster colors
library(ggplot2)
library(paletteer)
sc1000_LAD_Cryo_Sham_nonBlymphocytes1_fcol_original <- sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol <- as.character(paletteer_c("ggthemes::Green-Blue Diverging", 18))
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[17] <- "#79D4B8"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[15] <- "#7663A2"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[14] <- "#7663A2"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[11] <- "#FF8B00"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[18] <- "#FF0C00"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[1] <- "#0084A7"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[10] <- "#D083AF"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[6] <- "#C26395"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[4] <- "#850F27"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[16] <- "#E3AFDF"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[13] <- "#EAA928"

plotmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, um=T, cex=0.5)




# select and assign new cluster, medoids and color for gdT cell-like CD3e+ cells (CD4- CD8- CD44+):
require(Seurat)
library(patchwork)
library(ggplot2)
plot<-ggplot(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@umap, aes(V1, V2)) +
  geom_point()
select.cells <- CellSelector(plot = plot)
select.cells2 <- CellSelector(plot = plot)
select.cells3 <- CellSelector(plot = plot)
select.cells <- unique(c(select.cells, select.cells2, select.cells3))

cpart <- sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart
cpart[select.cells] <- max(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart) + 1
unique(cpart)
head(cpart)
head(cpart[select.cells])
m <- compmedoids(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,cpart)
fcol <- sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol
fcol[length(fcol)+1] <- "#26456E"
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart <- cpart
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol <- fcol
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@medoids <- m
plotmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, um=T, cex = 0.5)
table(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)








require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata), 7,8)
head(types)
unique(types)
types[types%in%"LA"] <- "LAD"
types[types%in%"MI"] <- "Cryoablation"
types[types%in%"Sh"] <- "Sham"
plotsymbolsmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, types, um=T, cex=0.5, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, types, um=T, subset = "LAD")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, types, um=T, subset = "Cryoablation")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, types, um=T, subset = "Sham", cex=1, samples_col  = rep("#0A1722", length(types[types%in%"Sham"])))

# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata), 7,12)

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
plotsymbolsmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, types, um=T, cex = 0.5, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))




plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Cd3e",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Cd4",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Cd8a",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Cd44",logsc=T,fr=F, um=T, cex=1) # T cell activation marker
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Sell",logsc=T,fr=F, um=T, cex=1) # CD62L, naïve T cells
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Foxp3",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Trdc",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Gzma",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Il2ra",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Il17a",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Mki67",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Trac",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Xcl1",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Pxdc1",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Il5",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Ifngr2",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Camk1d",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Cd14",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Isg15",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Klra5",logsc=T,fr=F, um=T, cex=1)

dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,c(11,18),pvalue=0.01)








# plot umap with cluster names annotated

cluster_ann <- sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart
cluster_ann[cluster_ann%in%1] <- "ILC2"
cluster_ann[cluster_ann%in%2] <- "T_CD4_effector"
cluster_ann[cluster_ann%in%3] <- "T_CD4_naive" 
cluster_ann[cluster_ann%in%4] <- "NK_Gzma" 
cluster_ann[cluster_ann%in%5] <- "T_CD8_naive" 
cluster_ann[cluster_ann%in%6] <- "NK_Klra5"
cluster_ann[cluster_ann%in%7] <- "T_CD8_effector" 
cluster_ann[cluster_ann%in%8] <- "T_IFNg_naive"
cluster_ann[cluster_ann%in%9] <- "T_CD4_naive"
cluster_ann[cluster_ann%in%10] <- "T_gd"
cluster_ann[cluster_ann%in%11] <- "MAIT" 
cluster_ann[cluster_ann%in%12] <- "T_Macro"
cluster_ann[cluster_ann%in%13] <- "T_Isg15" 
cluster_ann[cluster_ann%in%14] <- "Treg"
cluster_ann[cluster_ann%in%15] <- "Treg"
cluster_ann[cluster_ann%in%16] <- "NK_T" 
cluster_ann[cluster_ann%in%17] <- "ILC2_IL5" 
cluster_ann[cluster_ann%in%18] <- "MAIT_IL17" 


cluster_ann_nonBlymphocytes <- cluster_ann
plotsymbolsmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, cluster_ann_nonBlymphocytes, um=T, cex=0.5, samples_col =sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol
)



nonBlymphocytes_populations <- cluster_ann_nonBlymphocytes

# matching cell type names with colors
sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol
head(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)

nonBlymphocytes_fcol <- sc1000_LAD_Cryo_Sham_nonBlymphocytes1@fcol[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart]
names(nonBlymphocytes_fcol) <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)
head(nonBlymphocytes_fcol)

head(nonBlymphocytes_populations)
identical(names(nonBlymphocytes_populations), names(nonBlymphocytes_fcol))

nonBlymphocytes_clusters_fcol <- data.frame(clusters=nonBlymphocytes_populations, colors=nonBlymphocytes_fcol)

saveRDS(nonBlymphocytes_clusters_fcol, "df_cluster_colors_lymphocytes.rds")








# plot cell cycle scores

S_score   <- colMeans(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)),])

# plot cell cycle scores
plotfeatmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, G2M_score, "G2M score", logsc=T, um=T, cex=0.5)
plotfeatmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, S_score, "S phase score", logsc=T, um=T, cex=0.5)



# Plot box plot of G2M genes across time points (among ILC2)
ILC2_cells <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(1,17)]
df <- colSums(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@counts)
sample1 <- intersect(sample1, ILC2_cells)
sample3 <- intersect(sample3, ILC2_cells)
sample7 <- intersect(sample7, ILC2_cells)
sample28 <- intersect(sample28, ILC2_cells)
sample56 <- intersect(sample56, ILC2_cells)
sampleSh <- intersect(sampleSh, ILC2_cells)
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
p<-ggplot(df, aes(x=Days_post_MI, y=log2(G2M), fill=Days_post_MI)) + geom_boxplot()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))


# Plot box plot of S phase genes across time points (among ILC2)
ILC2_cells <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(1,17)]
df <- colSums(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[cc_genes$s,]) * median(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@counts)
# pull out cells of different time points (in non-B lymphocyte compartment)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
# pull out MAIT/gdT cells from each time point
sample1 <- intersect(sample1, ILC2_cells)
sample3 <- intersect(sample3, ILC2_cells)
sample7 <- intersect(sample7, ILC2_cells)
sample28 <- intersect(sample28, ILC2_cells)
sample56 <- intersect(sample56, ILC2_cells)
sampleSh <- intersect(sampleSh, ILC2_cells)
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
# log scale
p<-ggplot(df, aes(x=Days_post_MI, y=log2(S), fill=Days_post_MI)) + geom_boxplot()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))

# linear scale
p<-ggplot(df, aes(x=Days_post_MI, y=S, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 10) +
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))



# Plot box plot of S phase genes across time points (among gdT/MAIT-like cells)

MAIT_gdT <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(10,11,18)]
df <- colSums(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[cc_genes$s,]) * median(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@counts)
# pull out cells of different time points (in non-B lymphocyte compartment)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
# pull out MAIT/gdT cells from each time point
sample1 <- intersect(sample1, MAIT_gdT)
sample3 <- intersect(sample3, MAIT_gdT)
sample7 <- intersect(sample7, MAIT_gdT)
sample28 <- intersect(sample28, MAIT_gdT)
sample56 <- intersect(sample56, MAIT_gdT)
sampleSh <- intersect(sampleSh, MAIT_gdT)
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
# log scale
p<-ggplot(df, aes(x=Days_post_MI, y=log2(S), fill=Days_post_MI)) + geom_boxplot()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))

# linear scale
p<-ggplot(df, aes(x=Days_post_MI, y=S, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 10) +
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))





# Plot box plot of S phase genes across time points (among gdT cells only)

gdTcells <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(10)]
df <- colSums(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[cc_genes$s,]) * median(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@counts)
# pull out cells of different time points (in non-B lymphocyte compartment)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
# pull out MAIT/gdT cells from each time point
sample1 <- intersect(sample1, gdTcells)
sample3 <- intersect(sample3, gdTcells)
sample7 <- intersect(sample7, gdTcells)
sample28 <- intersect(sample28, gdTcells)
sample56 <- intersect(sample56, gdTcells)
sampleSh <- intersect(sampleSh, gdTcells)
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
# log scale
p<-ggplot(df, aes(x=Days_post_MI, y=log2(S), fill=Days_post_MI)) + geom_boxplot()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))

# linear scale
p<-ggplot(df, aes(x=Days_post_MI, y=S, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 10) +
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))





# Plot box plot of S phase genes across time points (among MAIT cells only)

MAIT_cells <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(11,18)]
df <- colSums(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata[cc_genes$s,]) * median(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@counts)
# pull out cells of different time points (in non-B lymphocyte compartment)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
# pull out MAIT/gdT cells from each time point
sample1 <- intersect(sample1, MAIT_cells)
sample3 <- intersect(sample3, MAIT_cells)
sample7 <- intersect(sample7, MAIT_cells)
sample28 <- intersect(sample28, MAIT_cells)
sample56 <- intersect(sample56, MAIT_cells)
sampleSh <- intersect(sampleSh, MAIT_cells)
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
# log scale
p<-ggplot(df, aes(x=Days_post_MI, y=log2(S), fill=Days_post_MI)) + geom_boxplot()+
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))

# linear scale
p<-ggplot(df, aes(x=Days_post_MI, y=S, fill=Days_post_MI)) + geom_boxplot()+ ylim(0, 10) +
  stat_compare_means(method = "anova", label.y = 7.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all
p + scale_fill_manual(values=c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))





















# interested genes

nonB_Lymphocytes_populations <- cluster_ann_nonBlymphocytes
table(nonB_Lymphocytes_populations)
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,c(13),pvalue=0.01)
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Ccl2",logsc=T,fr=F, um=T)

genes <- c(
  "Cd3e", "Trac",
  "Sell", "Ifngr2", "Ccr7", "Dusp10",
  "Ifng", "Cd44", "Cd4", "Cd8a", "Isg15","Ifit1", "Ifit3", "Cxcl10",
  "Foxp3", "Il2ra", "Ctla4", "Itgav", "Lamc1",
  "Il23r", "Trdc", "Tcrg-C1", "Cd163l1", "Kcnk1", "Pxdc1", "Rorc", "Mmp25", "Il17a",
  "Camk1d", "Xcl1", "Klra5", "Gzma", "Klra8", "Klrg1",
  "Gata3", "Kit", "Il7r", "Dach2", "Bmp2", "Il1rl1", "Il5", "Calca", "Il13"
  
  
)

fractDotPlot(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, genes, samples = nonB_Lymphocytes_populations, subset=c(
  "T_IFNg_naive", "T_CD4_naive", "T_CD8_naive",
  "T_CD4_effector", "T_CD8_effector", "T_Isg15",
  "Treg","T_gd", "MAIT", "MAIT_IL17",
  "NK_T", "NK_Klra5", "NK_Gzma",
  "ILC2", "ILC2_IL5"
), zsc = T, cap = 2)














# Dynamics of gene expression
# gene expression in a certain condition only
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
# pull out cells of different time points (in non-B lymphocyte compartment)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
#sample <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D028|NM_Ad_LA_D56|NM_Ad_LA_D28|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
plotexpmap(sc1000_LAD_Cryo_Sham_nonBlymphocytes1,"Tnfsf11",logsc=T,fr=F, um=T, cex=1, cells = sample7)


# subset ILC2 and expression of gene of interest
df <- sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata["Tnfsf11",names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(1,17)])]
# summed gene expression (normalized) in each time point
# normalized against datasize (total no. of cells in the corresponding timepoint)

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
timepoints <- c("Sham", "D1",   "D3",   "D7",   "D28",  "D56" )
names(a)<-timepoints
barplot(a, xlab = "Days post-MI", ylab = "normalized gene expression",main = "Tnfsf11 expression in ILC2", col = c("#0A1722", "#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2"))



# general gdT cell population dynamics (% of all immune cells)
sample1 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@ndata))]

gdT_cells <- names(sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart)[sc1000_LAD_Cryo_Sham_nonBlymphocytes1@cpart%in%c(10)]
v1 <- length(intersect(sample1, gdT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v3 <- length(intersect(sample3, gdT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v7 <- length(intersect(sample7, gdT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v28 <- length(intersect(sample28, gdT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
v56 <- length(intersect(sample56, gdT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
vSh <- length(intersect(sampleSh, gdT_cells)) / length(colnames(sc1000_LAD_Cryo_Sham_Immune@ndata)[grep("^NM_Ad_Sh|NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Immune@ndata))]) *100
a<-c(vSh,v1,v3,v7,v28,v56)
names(a)<-timepoints
barplot(a, xlab = "Days post-MI", ylab = "percentage of gdT cells",main = "% of gdT cells in immune compartment", col = c("#CC33FF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000CC"))

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

