

# subcluster vascular cells (ECs, pericytes and SMCs)

Vascular <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)[sc1000_LAD_Cryo_Sham_harmony@cpart%in%c(23,24,14,12,15)]


sc1000_LAD_Cryo_Sham_Vascular <-SCseq(sc1000_LAD_Cryo_Sham_harmony@expdata[,Vascular])
sc1000_LAD_Cryo_Sham_Vascular<-filterdata(sc1000_LAD_Cryo_Sham_Vascular,mintotal=1000, FGenes=rownames(sc1000_LAD_Cryo_Sham_Vascular@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(sc1000_LAD_Cryo_Sham_Vascular@expdata))])
expData  <- getExpData(sc1000_LAD_Cryo_Sham_Vascular)
# with batch removal (harmony) according to sample batches
batches <- Vascular
names(batches) <- Vascular
batch1 <- batches[grep("D01|D03|D07|D28|MI_D56",batch1)]
batch2 <- batches[grep("LA_D56|Sh_D56",batches)]
batch1 <- replace(batch1,,"b1")
batch2 <- replace(batch2,,"b2")
batchesFB <- c(batch1, batch2)

# regress out influences of cell cycle genes
S_score   <- colMeans(sc1000_LAD_Cryo_Sham_Vascular@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_Vascular@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_Vascular@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_Vascular@ndata)),])
regVar <- data.frame(S_score=S_score, G2M_score=G2M_score)
rownames(regVar) <- colnames(expData)


res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesFB, bmethod = "harmony", regVar = regVar)
# no batch effect removal
#res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham_Vascular@expdata)[sc1000_LAD_Cryo_Sham_Vascular@genes,colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)]
noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=10)
sc1000_LAD_Cryo_Sham_Vascular <- updateSC(sc1000_LAD_Cryo_Sham_Vascular,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham_Vascular <- comptsne(sc1000_LAD_Cryo_Sham_Vascular)
sc1000_LAD_Cryo_Sham_Vascular <- compumap(sc1000_LAD_Cryo_Sham_Vascular, n_neighbors = 30, min_dist = 0.8)
plotmap(sc1000_LAD_Cryo_Sham_Vascular, um=T, cex=0.2)
plotmap(sc1000_LAD_Cryo_Sham_Vascular, um=F, cex = 0.2)
table(sc1000_LAD_Cryo_Sham_Vascular@cpart)

# visualize different conditions
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
typesC <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata), 7,8)
head(typesC)
unique(typesC)
typesC[typesC%in%"LA"] <- "LAD"
typesC[typesC%in%"MI"] <- "Cryoablation"
typesC[typesC%in%"Sh"] <- "Sham"
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesC, um=T, cex=0.2, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesC, um=T, subset = "LAD")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesC, um=T, subset = "Cryoablation")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesC, um=T, subset = "Sham", samples_col  = rep("#0A1722", length(typesC[typesC%in%"Sham"])))

# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
typesT <- str_sub(colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata), 7,12)
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
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesT, cex=0.5, um=T, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesT, um=T, cex=0.5, subset = "D01", samples_col  = rep("#FF2E00", length(typesT[typesT%in%"D01"])))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesT, um=T, cex=0.5, subset = "D03", samples_col  = rep("#FFB900", length(typesT[typesT%in%"D03"])))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesT, um=T, cex=0.5, subset = "D07", samples_col  = rep("#0E7535", length(typesT[typesT%in%"D07"])))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesT, um=T, cex=0.5, subset = "D28")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, typesT, um=T, cex=0.5, subset = "D56")


# plot cell cycling genes
S_score   <- colMeans(sc1000_LAD_Cryo_Sham_Vascular@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_Vascular@ndata)),])
G2M_score <- colMeans(sc1000_LAD_Cryo_Sham_Vascular@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_Vascular@ndata)),])
plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=0.6, n = "G2M score")

sample1 <- colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)[grep("^NM_Ad_MI_D01|NM_Ad_LA_D01",colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata))]
sample3 <- colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)[grep("^NM_Ad_MI_D03|NM_Ad_LA_D03",colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata))]
sample7 <- colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)[grep("^NM_Ad_MI_D07|NM_Ad_LA_D07",colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata))]
sample28 <- colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)[grep("^NM_Ad_MI_D28|NM_Ad_LA_D28",colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata))]
sample56 <- colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)[grep("^NM_Ad_MI_D56|NM_Ad_LA_D56",colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata))]
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata)[grep("^NM_Ad_Sh",colnames(sc1000_LAD_Cryo_Sham_Vascular@ndata))]

plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=1, cells = sampleSh, n = "G2M score Sham")
plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=1, cells = sample1, n = "G2M score Day 1")
plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=1, cells = sample3, n = "G2M score Day 3")
plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=1, cells = sample7, n = "G2M score Day 7")
plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=1, cells = sample28, n = "G2M score Day 28")
plotfeatmap(sc1000_LAD_Cryo_Sham_Vascular, G2M_score, logsc=T, um=T, cex=1, cells = sample56, n = "G2M score Day 56")

# Violin plot G2M genes across clusters
violinMarkerPlot(cc_genes$g2m,sc1000_LAD_Cryo_Sham_Vascular,set=c(8,1,3,4,16,5), ti = "G2M genes")

# gene expressions box plots

df_G2M <- data.frame(expression=colSums(sc1000_LAD_Cryo_Sham_Vascular@ndata[cc_genes$g2m,]* median(sc1000_LAD_Cryo_Sham_Vascular@counts)), populations = as.factor(Vascular_populations))
df_S <- data.frame(expression=colSums(sc1000_LAD_Cryo_Sham_Vascular@ndata[cc_genes$s,]* median(sc1000_LAD_Cryo_Sham_Vascular@counts)), populations = as.factor(Vascular_populations))

# box plot of G2M genes
# log-scale
p<-ggplot(df_G2M, aes(x=Vascular_populations, y=log10(expression), fill=populations)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),) +
  
  stat_compare_means(method = "anova", label.y = -2)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all

p


# linear scale
p<-ggplot(df_G2M, aes(x=Vascular_populations, y=expression, fill=populations)) + ylim(0,20) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),) +
  
  stat_compare_means(method = "anova", label.y = 18)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all

p

# box plot of S genes
# log-scale
p<-ggplot(df_S, aes(x=Vascular_populations, y=log10(expression), fill=populations)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p






# Plot violin plot of G2M genes across time points
Macrophages <- names(sc1000_LAD_Cryo_Sham_Vascular@cpart)[sc1000_LAD_Cryo_Sham_Vascular@cpart%in%c(8,7,6,13,2,1,20,22)]
df <- colSums(sc1000_LAD_Cryo_Sham_Vascular@ndata[cc_genes$g2m,]) * median(sc1000_LAD_Cryo_Sham_Vascular@counts)
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








# gene markers for EC populations:
# CapEC: Aqp1, Timd4, Aplnr, Lpl, Car4
# CapA-EC: Rbp7, Cxcl12, Btnl9, Hey1, Unc5b
# CapV-EC: Aqp1, Eln
# StrEC: Fos, Fosb, Hspa1a, Hspa1b, Egr1
# IntEC: Isg15, Ifit3, Ifit1, Ifit2, Rsad2
# AngEC: Ift122, Fscn1
# ArtEC: Fbln5, Clu, Sema3g, Fn1
# LymEC: Ccl21a, Mmrn1, Fgl2, Lcn2, Lyve1
# CyclEC: Hmgb2, Top2a, Cenpf, Stmn1, Hist1h2ap
# RepEC: Mcm5, Mcm3, Dut, Hells, Lig1
# EndoEC: Vwf, Vcam1, Mgp, Cfh, Apoe, Dcn, Igfbp5, Npr3

plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Cdh5",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Vwf",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Lyve1",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Aqp1",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Timp4",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Ptprc",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Rbp7",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Eln",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Fos",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Fosb",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Isg15",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Ift122",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Fscn1",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Fbln5",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Sema3g",logsc=T,fr=F, um=T, cex=1)



# gene markers for pericytes/SMC populations
# Peri1: Kcnj8, Vtn, Plat, Coro1b, Gas1
# Peri2: Mgp, Cxcl12, Ifitm1
# VSMC2: Cd36
# VSMC1: Acta2, Tagln, Myh11, Mustn1, Myl9
# StrVSMC: Ppp1r15a, Btg2, Rasd1, Nr4a1
# StrPeri: Hspa1a, Jun, Fosb
# IntPeri: Isg15, Ly6e, Ifit3, Ifit2, Iigp1

plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Pdgfrb",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Acta2",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Kcnj8",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Mgp",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Cxcl12",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Ifitm1",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Cd36",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Tagln",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Myh11",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Isg15",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Ly6e",logsc=T,fr=F, um=T, cex=1)
plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Ifit3",logsc=T,fr=F, um=T, cex=1)




A <- names(sc1000_LAD_Cryo_Sham_Vascular@cpart)[sc1000_LAD_Cryo_Sham_Vascular@cpart %in% c(5,3,2,6)]
B <- names(sc1000_LAD_Cryo_Sham_Vascular@cpart)[sc1000_LAD_Cryo_Sham_Vascular@cpart %in% c(7,9)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_Vascular,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="D56_Shared_EC",Bname="D56_LAD_EC",show_names=TRUE,padj=TRUE)

plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Cdh5",logsc=T,fr=F, um=T)

dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Vascular,c(19),pvalue=0.01)


# plot umap with cluster names annotated

cluster_ann <- sc1000_LAD_Cryo_Sham_Vascular@cpart
cluster_ann[cluster_ann%in%1] <- "vEC_capillary1"
cluster_ann[cluster_ann%in%2] <- "vSMC2"
cluster_ann[cluster_ann%in%3] <- "vEC_capillary2" 
cluster_ann[cluster_ann%in%4] <- "vEC_angio_IFN" 
cluster_ann[cluster_ann%in%5] <- "vEC_Endocardial" 
cluster_ann[cluster_ann%in%6] <- "vPericyte_INF"
cluster_ann[cluster_ann%in%7] <- "vPericyte_quiescent" 
cluster_ann[cluster_ann%in%8] <- "vEC_Arterial"
cluster_ann[cluster_ann%in%9] <- "vSMC2"
cluster_ann[cluster_ann%in%10] <- "vEC_Lymphatic"
cluster_ann[cluster_ann%in%11] <- "vSMC1" 
cluster_ann[cluster_ann%in%12] <- "vEC_Immune"
cluster_ann[cluster_ann%in%13] <- "vPericyte_FB" 
cluster_ann[cluster_ann%in%14] <- "vSMC_Ryr2"
cluster_ann[cluster_ann%in%15] <- "vEC_Immune"
cluster_ann[cluster_ann%in%16] <- "vEC_Areg_Dkk2_Wnt" 
cluster_ann[cluster_ann%in%17] <- "vEC_metabolic" 
cluster_ann[cluster_ann%in%18] <- "vEC_Immune" 
cluster_ann[cluster_ann%in%19] <- "vEpicardial_derived" 

cluster_ann_Vascular <- cluster_ann

plotsymbolsmap(sc1000_LAD_Cryo_Sham_Vascular, cluster_ann_Vascular, um=T, cex=0.2, samples_col = 
                 sc1000_LAD_Cryo_Sham_Vascular@fcol
)

Vascular_populations <- cluster_ann_Vascular

# matching cell type names with colors
sc1000_LAD_Cryo_Sham_Vascular@fcol
head(sc1000_LAD_Cryo_Sham_Vascular@cpart)

Vascular_fcol <- sc1000_LAD_Cryo_Sham_Vascular@fcol[sc1000_LAD_Cryo_Sham_Vascular@cpart]
names(Vascular_fcol) <- names(sc1000_LAD_Cryo_Sham_Vascular@cpart)
head(Vascular_fcol)

head(Vascular_populations)
identical(names(Vascular_populations), names(Vascular_fcol))

Vascular_clusters_fcol <- data.frame(clusters=Vascular_populations, colors=Vascular_fcol)

saveRDS(Vascular_clusters_fcol, "df_cluster_colors_Vascular.rds")








# dotplots of selected genes
# pool clusters into defined sub-populations
plotmap(sc1000_LAD_Cryo_Sham_Vascular, um=T, cex=0.5)
Vascular_populations <- cluster_ann_Vascular
table(Vascular_populations)

# interested genes
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_Vascular,c(2,3),pvalue=0.01)

genes <- c("Cdh5", "Fbln5", "Sema3g", # Aterial EC
           "Rbp7", "Cxcl12", "Aqp1", "Timp4", "Isg15", "Ifit3", "Ift122", "Fscn1","Prnp", "Areg", "Dkk2", # Capillary EC
           "Vwf", "Vcam1", # Endocardial EC
           "Lyve1", "Mmrn1", # Lymphatic EC
           "Hexb", "Camk1d", # Metabolic EC
           "Msln", # Epicardial-derived
           "Pdgfrb", "Mgp", "Kcnj8", "Pdgfra", # pericytes
           "Acta2","Myh11", "Tcap", "Slc38a11","Ryr2" # SMC
           
)

fractDotPlot(sc1000_LAD_Cryo_Sham_Vascular, genes, samples = Vascular_populations, subset=c(
  "vEC_Arterial","vEC_capillary1", "vEC_capillary2", "vEC_angio_IFN", "vEC_Areg_Dkk2_Wnt",
  "vEC_Endocardial", "vEC_Lymphatic", "vEC_metabolic",  
  "vEpicardial_derived", 
  "vPericyte_quiescent","vPericyte_FB", "vPericyte_INF",
  "vSMC1", "vSMC2","vSMC_Ryr2"
), zsc = T, cap = 2)


plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Ccl2",logsc=T,fr=F, um=T)



plotexpmap(sc1000_LAD_Cryo_Sham_Vascular,"Axl",logsc=T,fr=F, um=T)












# Transition probabilities
expData  <- getExpData(sc1000_LAD_Cryo_Sham_Vascular)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=6,seed=12345, batch = batchesFB, bmethod = "harmony", regVar = regVar)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
# update cl after assigning new cluster identities
cl$partition <- sc1000_LAD_Cryo_Sham_Vascular@cpart

probs_Vas <-transitionProbs(res,cl,pvalue=0.01)
plotTrProbs(sc1000_LAD_Cryo_Sham_Vascular,probs_Vas,tp=.5,prthr=0,cthr=0,fr=F, um = T, cex=0.1)

# select clusters for pseudo-temporal plot (Slingshot)
require(FateID)
require(SingleCellExperiment)
require(slingshot)
require(DelayedMatrixStats)
set <- c(8,3,1,4) # for EC angiogenesis lineage
set <- c(7,6,19,13) # for pericyte/Epicardial/FB lineage
set <- c(7,6,11,2) # for Pericyte-SMC lineage
set <- c(11,2,14) # for SMC-CM-like lineage

pt <- pseudoTime(sc1000_LAD_Cryo_Sham_Vascular,m="umap",set=set)
plotPT(pt,sc1000_LAD_Cryo_Sham_Vascular,clusters=FALSE)
plotPT(pt,sc1000_LAD_Cryo_Sham_Vascular,clusters=T)

# plot gene module heat map
fs <- extractCounts(sc1000_LAD_Cryo_Sham_Vascular,minexpr=5,minnumber=5,pt=pt)
s1d   <- getsom(fs,nb=50,alpha=1)
ps    <- procsom(s1d,corthr=.85,minsom=0)
part  <- pt$part
ord   <- pt$ord
plotheatmap(ps$all.z, xpart=part[ord], xcol=sc1000_LAD_Cryo_Sham_Vascular@fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)
# gene list in each module:
module_genelist <- ps$nodes
names(module_genelist)[module_genelist%in%c(9)] # input module number
modules_dediff <- names(module_genelist)[module_genelist%in%c(9,12,14,15)]
# plotting only gene modules 1-9:
plotheatmap(ps$all.z[names(module_genelist)[module_genelist%in%c(9:13)],], xpart=part[ord], xcol=sc1000_LAD_Cryo_Sham_Vascular@fcol, ypart=ps$nodes[module_genelist%in%c(9:13)], xgrid=FALSE, ygrid=TRUE, xlab=TRUE)

# plot selected gene expression
plotexpression(fs,y=part,g="Erc2",n=ord,col=sc1000_LAD_Cryo_Sham_Vascular@fcol,cex=1,alpha=1)














