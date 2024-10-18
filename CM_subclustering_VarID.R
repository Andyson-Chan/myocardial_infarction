# subclustering CM



# prepare for Harmony batch effect removal (perform in "res"), according to sample harvesting batches
require(stringr)
batch1 <- rep("b1", dim(prdata_LAD_CM)[2])
names(batch1) <- colnames(prdata_LAD_CM)

batch2 <- rep("b2", dim(cbind(prdata_Cr_CM_Day1_1, prdata_Cr_CM_Day3_1, prdata_Cr_CM_Day7_1, prdata_Cr_CM_Day28_1, prdata_Cr_CM_Day56_1))[2])
names(batch2) <- colnames(cbind(prdata_Cr_CM_Day1_1, prdata_Cr_CM_Day3_1, prdata_Cr_CM_Day7_1, prdata_Cr_CM_Day28_1, prdata_Cr_CM_Day56_1))

batch3 <- rep("b3", dim(cbind(prdata_Cr_CM_Day1_2, prdata_Cr_CM_Day3_2, prdata_Cr_CM_Day7_2, prdata_Cr_CM_Day28_2, prdata_Cr_CM_Day56_2))[2])
names(batch3) <- colnames(cbind(prdata_Cr_CM_Day1_2, prdata_Cr_CM_Day3_2, prdata_Cr_CM_Day7_2, prdata_Cr_CM_Day28_2, prdata_Cr_CM_Day56_2))

batch4 <- rep("b4", dim(cbind(prdata_Sham_CM_Day1, prdata_Sham_CM_Day3, prdata_Sham_CM_Day7, prdata_Sham_CM_Day28, prdata_Sham_CM_Day56))[2])
names(batch4) <- colnames(cbind(prdata_Sham_CM_Day1, prdata_Sham_CM_Day3, prdata_Sham_CM_Day7, prdata_Sham_CM_Day28, prdata_Sham_CM_Day56))



batchesCM_LAD_Cryo <- c(batch1, batch2, batch3, batch4)

length(batchesCM_LAD_Cryo)

length(unique(names(batchesCM_LAD_Cryo)))

# perform clustering with VarID (decided to use VarDI batch effect removal)

scCM_LAD_Cryo_sham1 <-SCseq(prdata_CM_LAD_Cryo_Sham)
scCM_LAD_Cryo_sham1<-filterdata(scCM_LAD_Cryo_sham1,mintotal=1000, FGenes=rownames(scCM_LAD_Cryo_sham1@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(scCM_LAD_Cryo_sham1@expdata))])
# remove also "Cmss1" which is up in both LAD and Sham datasets (technical)
#scCM_LAD_Cryo_sham1<-filterdata(scCM_LAD_Cryo_sham1,mintotal=1000, CGenes=rownames(scCM_LAD_Cryo_sham1@expdata)[grep("^(mt|Rp(l|s)|Gm\\d|Cmss1)",rownames(scCM_LAD_Cryo_sham1@expdata))])
expData  <- getExpData(scCM_LAD_Cryo_sham1)
res   <- pruneKnn(expData, knn=7, no_cores=8,batch=batchesCM_LAD_Cryo)
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
scCM_LAD_Cryo_sham1 <- updateSC(scCM_LAD_Cryo_sham1,res=res,cl=cl,flo=.1)
scCM_LAD_Cryo_sham1 <- comptsne(scCM_LAD_Cryo_sham1)
scCM_LAD_Cryo_sham1 <- compumap(scCM_LAD_Cryo_sham1, spread = 7, min_dist = 0.5)
plotmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5)
plotmap(scCM_LAD_Cryo_sham1, um=F, cex=0.5)
table(scCM_LAD_Cryo_sham1@cpart)

# visualize cells of different time points
# examples of names: CM1Ad_Cr_D56, CM1Ad_LA_D28, CM1Ad_Sh_D28
types <- str_sub(colnames(scCM_LAD_Cryo_sham1@ndata), 7,12)
# merge all Sham time points into 1 type
types<- gsub("Sh_D01", "Sham", types)
types<- gsub("Sh_D03", "Sham", types)
types<- gsub("Sh_D07", "Sham", types)
types<- gsub("Sh_D28", "Sham", types)
types<- gsub("Sh_D56", "Sham", types)

types<- gsub("LA_D01", "D01", types)
types<- gsub("LA_D03", "D03", types)
types<- gsub("LA_D07", "D07", types)
types<- gsub("LA_D28", "D28", types)
types<- gsub("LA_D56", "D56", types)

types<- gsub("Cr_D01", "D01", types)
types<- gsub("Cr_D03", "D03", types)
types<- gsub("Cr_D07", "D07", types)
types<- gsub("Cr_D28", "D28", types)
types<- gsub("Cr_D56", "D56", types)
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, cex = 0.5, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "D01")
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "D03")
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "D07")
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "D28")
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "D56")
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "Sham")

# visualize cells of different conditions
# examples of names: CM1Ad_Cr_D56, CM1Ad_MI_D28, NM_Ad_Sh_D01
types <- str_sub(colnames(scCM_LAD_Cryo_sham1@ndata), 7,8)
types<- gsub("Cr", "Cryoablation", types)
types<- gsub("LA", "LAD", types)
types<- gsub("Sh", "Sham", types)
head(types)
unique(types)
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, cex=0.5, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, subset = "Sham", cex=1, samples_col  = rep("#0A1722", length(types[types%in%"Sham"])))

# visualize cells of different time points and conditions
# examples of names: CM1Ad_Cr_D56, CM1Ad_LA_D28, CM1Ad_Sh_D28
require(stringr)
types <- str_sub(colnames(scCM_LAD_Cryo_sham1@ndata), 7,12)
head(types)
unique(types)
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, cex=1)
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, cex=1, subset = c("LA_D01", "LA_D03", "LA_D07", "LA_D28", "LA_D56"))
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, cex=1, subset = c("Cr_D01", "Cr_D03", "Cr_D07", "Cr_D28", "Cr_D56"))
plotsymbolsmap(scCM_LAD_Cryo_sham1, types, um=T, cex=1, subset = c("Sh_D01", "Sh_D03", "Sh_D07", "Sh_D28", "Sh_D56"))



plotexpmap(scCM_LAD_Cryo_sham1,"Nr4a1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Mb",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Nppa",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Angpt1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Cmss1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Slit2",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Ppargc1a",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Ddx60",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Btk",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Ankrd1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Gpc6",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Rasef",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Gck",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Fbn2",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Dlgap1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Fgf1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Hif1a",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Vegfb",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Myh6",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Myh7",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,"Yap1",logsc=T,fr=F, um=T, cex=0.5)
plotexpmap(scCM_LAD_Cryo_sham1,c("Myh7", "Actc1", "Mb", "Mdh2"), n = "Dediff genes",logsc=T,fr=F, um=T, cex=1)


plotfeatmap(scCM_LAD_Cryo_sham1, scCM_LAD_Cryo_sham1@expdata["Mb",], n = "UMI Mb",logsc=F, cex=1, ceil = 40)
plotfeatmap(scCM_LAD_Cryo_sham1, scCM_LAD_Cryo_sham1@expdata["Ckm",], n = "UMI Ckm",logsc=F, cex=1, ceil = 10)
plotfeatmap(scCM_LAD_Cryo_sham1, scCM_LAD_Cryo_sham1@expdata["Eno3",], n = "UMI Eno3",logsc=F, cex=1, ceil = 6)



# plot cell cycle scores
S_score   <- colMeans(scCM_LAD_Cryo_sham1@ndata[intersect(cc_genes$s,rownames(sc1000_LAD_Cryo_Sham_Immune1@ndata)),])
G2M_score <- colMeans(scCM_LAD_Cryo_sham1@ndata[intersect(cc_genes$g2m,rownames(sc1000_LAD_Cryo_Sham_Immune1@ndata)),])

plotfeatmap(scCM_LAD_Cryo_sham1, G2M_score, "G2M score", logsc=F, um=T, cex=2)
plotfeatmap(scCM_LAD_Cryo_sham1, S_score, "S phase score", logsc=F, um=T, cex=2)
plotexpmap(scCM_LAD_Cryo_sham1, cc_genes$s, "S phase score", um=T, logsc=T, cex=2)
plotexpmap(scCM_LAD_Cryo_sham1, cc_genes$g2m, "G2M phase score", um=T, logsc=T, cex=2)


violinMarkerPlot(cc_genes$g2m,scCM_LAD_Cryo_sham1, ti = "G2M genes",set=c(5,2,10,12,17,15,14))
violinMarkerPlot(cc_genes$s,scCM_LAD_Cryo_sham1, ti = "S genes",set=c(5,2,10,12,17,15,14))

violinMarkerPlot(c("Pcna", "Mcm2"),scCM_LAD_Cryo_sham1, ti = "Pcna and Mcm2",set=c(5,2,10,12,17,15,14))






# collect CM cell identities of different populations
CM_sham <- colnames(scCM_LAD_Cryo_sham1@ndata)[grep("Sh",colnames(scCM_LAD_Cryo_sham1@ndata))]
CM_normal <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(4,5:9,11,13,16)])
CM_prehyper <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(1,12)])
CM_hyper <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(17)])
CM_dediff <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(15)])
CM_angio <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(2,3)])
CM_slit2 <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(10)])
CM_Btk <- names(scCM_LAD_Cryo_sham1@cpart[scCM_LAD_Cryo_sham1@cpart%in%c(14)])



# dotplots of selected genes
# pool clusters into defined sub-populations
plotmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5)
#CM_populations <- as.character(scCM_LAD_Cryo_sham1@cpart)
CM_populations <- scCM_LAD_Cryo_sham1@cpart
CM_populations[CM_populations %in% c(1,4:9,11,13,16)] <- "Homeostatic"
CM_populations[CM_populations %in% c(14)] <- "Btk+"
CM_populations[CM_populations %in% c(12)] <- "Pre-hypertrophic"
CM_populations[CM_populations %in% c(17)] <- "Hypertrophic"
CM_populations[CM_populations %in% c(15)] <- "Dedifferentiating"
CM_populations[CM_populations %in% c(2,3)] <- "Angiogenic"
CM_populations[CM_populations %in% c(10)] <- "Slit2+"

table(CM_populations)
# interested genes
dg <- clustdiffgenes(scCM_LAD_Cryo_sham1,c(2,3),pvalue=0.01)

genes = c("Myh7", "Ankrd1", "Mybpc2", "Acta1", "Osmr",   # pre-hypertrophic
          "Nppa", "Nppb","Gdf15", "Nrxn3", "Xirp2",      # hypertrophic
          "Actc1","Mb", "Mdh2", "Cox6a2", "Atp5e",       # dedifferentiating CM
          "Angpt1", "Slit2", "Slit3","Fgf12", "Cntn2",  # angiogenic
          "Btk" # Btk+
)
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Acta1")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Xirp2")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Nppb")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Myh7")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Btk")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Angpt1")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Ryr2")
plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex=0.5, logsc = T, "Il31ra")
fractDotPlot(scCM_LAD_Cryo_sham1, genes, samples = CM_populations, logscale=TRUE)
fractDotPlot(scCM_LAD_Cryo_sham1, genes, samples = CM_populations, subset=c("Homeostatic","Btk+","Slit2+","Angiogenic","Dedifferentiating","Hypertrophic","Pre-hypertrophic"), zsc = T, cap = 1.5)
table(CM_populations)

# dedifferentiation-related receptor genes
genes2 <- c("Erbb2", "Osmr", "Fgfr1", "Tnfrsf12a", "Notch1", "Fzd3", "Il1r1", "Bmpr2", "Ncl", "Sdc2", "Pecam1", "Cdh5", "Pdgfra", "Pdgfrb", "Itgb1","Kit", "Ly6a", "Kdr", "Gpc1", "Cd47", "Flt1", "Des", "Ptprj", "Bcam")
fractDotPlot(scCM_LAD_Cryo_sham1, genes2, samples = CM_populations, subset=c("Homeostatic","Btk+","Slit2+","Angiogenic","Dedifferentiating","Hypertrophic","Pre-hypertrophic"), zsc = T, cap = 1.5)



# gene expressions box plots

df_G2M <- data.frame(expression=colSums(scCM_LAD_Cryo_sham1@ndata[cc_genes$g2m,]), populations = as.factor(CM_populations))
df_S <- data.frame(expression=colSums(scCM_LAD_Cryo_sham1@ndata[cc_genes$s,]), populations = as.factor(CM_populations))

# box plot of G2M genes
# log-scale
p<-ggplot(df_G2M, aes(x=CM_populations, y=log10(expression), fill=populations)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),) +
  
  stat_compare_means(method = "anova", label.y = -2)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   # Pairwise comparison against all

p





# box plot of S genes
# log-scale
p<-ggplot(df_S, aes(x=CM_populations, y=log10(expression), fill=populations)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p








# heatmaps of genee expressions
genes <- c("Nkx2-5", "Actc1", "Mb", "Myh6", "Il1r1", "Bmpr2", "Tnfrsf12a")

fractDotPlot(scCM_LAD_Cryo_sham1, genes, cl=c(15, 11, 4, 8, 5, 9, 2, 13, 6, 10, 17), zsc=TRUE)

