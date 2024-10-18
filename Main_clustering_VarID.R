reticulate::use_python("/python", required = T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
modules
reticulate::py_available()
require(RaceID)
require(Matrix)
require(ggplot2)
require(stringr)

# load gene annotated gene expression matrices:
# prdata_LAD, prdata_Cryoablation, prdata_sham


# prepare for Harmony batch effect removal (perform in "res")
# LAD and Sham were performed in the same batch; Cryoablation is performed in another batch
batch1 <- rep("b1", dim(cbind(prdata_LAD, prdata_sham))[2])
names(batch1) <- colnames(cbind(prdata_LAD, prdata_sham))
batch2 <- rep("b2", dim(prdata_Cryoablation)[2])
names(batch2) <- colnames(prdata_Cryoablation)
batches <- c(batch1, batch2)
length(batches)

# perform clustering
sc1000_LAD_Cryo_Sham <-SCseq(prdata_LAD_Cryo_Sham)
sc1000_LAD_Cryo_Sham<-filterdata(sc1000_LAD_Cryo_Sham,mintotal=1000, FGenes=rownames(sc1000_LAD_Cryo_Sham@expdata)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(sc1000_LAD_Cryo_Sham@expdata))])
expData  <- getExpData(sc1000_LAD_Cryo_Sham)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,no_cores=8,seed=12345, FSelect = T,batch = batches, bmethod = "harmony")
cl    <- graphCluster(res,pvalue=0.01, use.weights = T, use.leiden = T, leiden.resolution = 1)
table(cl$partition)
probs <- transitionProbs(res,cl)
x     <- as.matrix(sc1000_LAD_Cryo_Sham@expdata)[sc1000_LAD_Cryo_Sham@genes,colnames(sc1000_LAD_Cryo_Sham@ndata)]
sc1000_LAD_Cryo_Sham <- updateSC(sc1000_LAD_Cryo_Sham,res=res,cl=cl,flo=.1)
sc1000_LAD_Cryo_Sham <- comptsne(sc1000_LAD_Cryo_Sham)
sc1000_LAD_Cryo_Sham <- compumap(sc1000_LAD_Cryo_Sham, spread = 1, min_dist = 0.5)
plotmap(sc1000_LAD_Cryo_Sham, um=T)
plotmap(sc1000_LAD_Cryo_Sham, um=F, cex = 0.2)
table(sc1000_LAD_Cryo_Sham@cpart)
plotmap(sc1000_LAD_Cryo_Sham_harmony, um=F, cex = 0.2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# visualize different conditions
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 7,8)
head(types)
unique(types)
types[types%in%"LA"] <- "LAD"
types[types%in%"MI"] <- "Cryoablation"
types[types%in%"Sh"] <- "Sham"
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, cex=0.2, samples_col = c("#008BCC","#C9655E","#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, cex=0.2, subset = "LAD")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, cex=0.2, subset = "Cryoablation")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, cex=0.2, subset = "Sham")
# visualize NM and CM
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 1,2)
head(types)
unique(types)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "CM", cex=0.2)
# visualize cells of different time points
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 7,12)
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

types<- gsub("MI_D01", "D01", types)
types<- gsub("MI_D03", "D03", types)
types<- gsub("MI_D07", "D07", types)
types<- gsub("MI_D28", "D28", types)
types<- gsub("MI_D56", "D56", types)


head(types)
unique(types)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, cex=0.2, samples_col = c("#FF2E00", "#FFB900", "#0E7535", "#0099FF", "#5B4DA2", "#0A1722"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "D01", cex=0.2)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "D03", cex=0.2)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "D07", cex=0.2)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "D28", cex=0.2)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "D56", cex=0.2)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, types, um=F, subset = "Sham", cex=0.2)

# create vectors with all cells and corresponding time points
cells_timepoints <- types
names(cells_timepoints) <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)
head(cells_timepoints)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cells_timepoints, um=F, cex=0.1)

# dotplots of marker gene expressions across major cell types
# first define clusters into major cell types
cluster_ann_major <- sc1000_LAD_Cryo_Sham_harmony@cpart
plotmap(sc1000_LAD_Cryo_Sham_harmony, um=F, cex=0.2)
# clusters:
# CM (4,27)
cluster_ann_major[cluster_ann_major%in%c(4,27)] <- "CM"
# quiescent FB (3,1)
cluster_ann_major[cluster_ann_major%in%c(3,1)] <- "quiescent FB"
# activated/myoFB (17,21,2,13,26)
cluster_ann_major[cluster_ann_major%in%c(17,21,2,13,26)] <- "activated/myoFB"
# Neutrophil (11,9)
cluster_ann_major[cluster_ann_major%in%c(11,9)] <- "Neutrophil"
# Macrophages (18,8,6,10,7,28,30)
cluster_ann_major[cluster_ann_major%in%c(18,8,6,10,7,28,30)] <- "Macrophages"
# DC (20,5)
cluster_ann_major[cluster_ann_major%in%c(20,5)] <- "DC"
# Mast_cell (31)
cluster_ann_major[cluster_ann_major%in%c(31)] <- "Mast_cells"
# EC_endocardial (24) 
cluster_ann_major[cluster_ann_major%in%c(24)] <- "EC_endocardial"
# EC_capillary (12,15) 
cluster_ann_major[cluster_ann_major%in%c(12,15)] <- "EC_capillary"
# Pericyte/SMC (14,23)
cluster_ann_major[cluster_ann_major%in%c(14,23)] <- "Pericyte/SMC"
# B_cells (22)
cluster_ann_major[cluster_ann_major%in%c(22)] <- "B_cells"
# T/NK_cell (16, 19) 
cluster_ann_major[cluster_ann_major%in%c(16,19)] <- "T/NK_cell"
# ILC2 (25)
cluster_ann_major[cluster_ann_major%in%c(25)] <- "ILC2"
# SwC (29)
cluster_ann_major[cluster_ann_major%in%c(29)] <- "SwC"


unique(cluster_ann_major)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_major, um=F, cex=0.1)

plotexpmap(sc1000_LAD_Cryo_Sham_harmony, um=F, logsc=T, cex=0.5, "S100a1")

# 14 major populations

genes <- c(
  "Myh6", "Tnnt2","Ttn","Gata6", #CM
  "Pdgfra", "Tcf21","Duox1", # quiescent FB
  "Postn", "Acta2","Ccl2", "Thbs4", # activated/myoFB
  "Ptprc", "H2-Ab1",# all immune
  "S100a8", "Ly6g", #Neutrophils
  "Ccr2", "Itgam","Cx3cr1","Timd4", #Macrophages
  "Itgax", "Cd209a", "Itgae", #DC
  "Cd3e","Cd4","Cd8a","Gzmb", #T/NK cells
  "Gata3","Il7r","Kit", #ILC2
  "Cd19", "Ly6d", #B cells
  "Cpa3", "Cma1", #Mast cells
  "Pdgfrb", "Myh11", "Kcnj8", #pericytes/SMC
  "Cdh5", "Pecam1", "Fbln5", "Cxcl12", #EC 
  "Vwf", "Vcam1", "Mgp", # EC endocardial
  "Plp1", "Prnp", "Ncam1", "S100b" # SwC
)

fractDotPlot(sc1000_LAD_Cryo_Sham_harmony, genes, samples = cluster_ann_major, zsc = F, subset = 
               c( "CM", 
                  "quiescent FB", "activated/myoFB",
                  "Neutrophil", "Macrophages", "DC",             
                  "T/NK_cell", "ILC2", "B_cells", "Mast_cells",
                  "Pericyte/SMC", "EC_capillary","EC_endocardial", "SwC" ) # align the cell type orders
)  # ,cap = 2




merscope_gene_markers<- c('Pax5', 
                          'Cd3e','Cd4','Gzmb','Tbx21', 
                          'Nppb','Alpl',
                          'Des','Ptgds','Nkx2-5','Tnfrsf12a',
                          'Mdh2','Ephx2','Hand2', 'Nr4a1',
                          'Ccr7', 'Mreg', 'Ccl22',
                          'Snx22','Xcr1','Klri1','Cd209a',
                          'Klk1','Siglech','Lag3', 'Dntt',
                          'Mcam','Aplnr','Fscn1','Nes','Kdr',
                          'Ckap2l','Mxd3','Gypa', 'Arg1',
                          'Sox9','Lox','Mfap4','Pdgfrb','Fbln2','Cd248',  'Axl', 'Sema3c','Sema3d',
                          'Ebf1',
                          'Gata3','Il17rb','Il2ra',
                          'Rorc','Kcnk1','Il18r1',
                          'Csf1r','Gpnmb','Fcgr3','Mrc1',   'Cx3cr1','Gas6','Pros1',
                          'Ifngr1','Hif1a','Ccr2','Cd14','Itgam',
                          'Cd163','Stab1','Tnfsf12'    ,'Timd4','Lyve1',
                          'Cd200r3', 'Kit', 'Il4',
                          'Tbx21','Xcl1','Klrb1c','Nkg7',
                          'Hp','Mmp9',
                          'Vtn','Ednra','Tbx3',
                          'Gja5','Igf2','Nrip2',
                          'Gfra3','Sox10',
                          'Tnfsf11',
                          'Foxp3','Ctla4',
                          'Flt4','Ndrg1','Adgrg3','Mmrn1',
                          'Msln','Krt8','Krt19','Krt7',
                          'Mki67','Pcna','Aurkb','Hmmr')

plotexpmap(scCM_LAD_Cryo_sham1, um=T, cex = 1, logsc=T, merscope_gene_markers)



# inspect UMI count and gene count in each sample

# plot normalized expressions across conditions
umi <- colSums(sc1000_LAD_Cryo_Sham_harmony@expdata)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@expdata), 7,12)
# merge all Sham time points into 1 type
types<- gsub("Sh_D01", "Sham D01", types)
types<- gsub("Sh_D03", "Sham D03", types)
types<- gsub("Sh_D07", "Sham D07", types)
types<- gsub("Sh_D28", "Sham D28", types)
types<- gsub("Sh_D56", "Sham D56", types)

types<- gsub("LA_D01", "LAD D01", types)
types<- gsub("LA_D03", "LAD D03", types)
types<- gsub("LA_D07", "LAD D07", types)
types<- gsub("LA_D28", "LAD D28", types)
types<- gsub("LA_D56", "LAD D56", types)

types<- gsub("MI_D01", "Cryo D01", types)
types<- gsub("MI_D03", "Cryo D03", types)
types<- gsub("MI_D07", "Cryo D07", types)
types<- gsub("MI_D28", "Cryo D28", types)
types<- gsub("MI_D56", "Cryo D56", types)

df_umi <- data.frame(umi, types)
df_umi <- df_umi[grep( "MI|LA|Sh_D01|Sh_D03|Sh_D28|Sh_D56" ,colnames(sc1000_LAD_Cryo_Sham_harmony@expdata)), ]
head(df_umi)
# violin plot
p<-ggplot(df_umi, aes(x=types, y=log10(umi), fill=types)) +
  geom_violin(trim=T)
p

# box plot
p<-ggplot(df_umi, aes(x=types, y=log10(umi), fill=types)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p

p<-ggplot(df_umi, aes(x=types, y=umi, fill=types)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p


# gene count
genecount <- colSums(sc1000_LAD_Cryo_Sham_harmony@expdata > 0)
df_gene <- data.frame(genecount, types)
df_gene <- df_gene[grep( "MI|LA|Sh_D01|Sh_D03|Sh_D28|Sh_D56" ,colnames(sc1000_LAD_Cryo_Sham_harmony@expdata)), ]

# box plot
p<-ggplot(df_gene, aes(x=types, y=log10(genecount), fill=types)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p
p<-ggplot(df_gene, aes(x=types, y=genecount, fill=types)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p


# mt %
mt_umi_count <- colSums(sc1000_LAD_Cryo_Sham_harmony@expdata[grep("^mt-", rownames(sc1000_LAD_Cryo_Sham_harmony@expdata)),])
mt_umi_count <- mt_umi_count / umi *100
df_mt <- data.frame(mt_umi_count, types)
df_mt <- df_mt[grep( "MI|LA|Sh_D01|Sh_D03|Sh_D28|Sh_D56" ,colnames(sc1000_LAD_Cryo_Sham_harmony@expdata)), ]

# box plot
# log-scale
p<-ggplot(df_mt, aes(x=types, y=log10(mt_umi_count), fill=types)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p
# linear scale
p<-ggplot(df_mt, aes(x=types, y=mt_umi_count, fill=types)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 8),)
p



# hypoxic gene expressions

genelist <- fgseaRes[fgseaRes$pathway%in%"Cellular response to hypoxia",leadingEdge] # from REpathway (see GESA script)
genelist <- bitr(genelist[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
genelist
Hypoxic_genes<- genelist$SYMBOL


Hypoxic_exp <- colSums(sc1000_LAD_Cryo_Sham_harmony@ndata[Hypoxic_genes,])


sampleLAD <- colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("LA_",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))]
plotfeatmap(sc1000_LAD_Cryo_Sham_harmony, Hypoxic_exp, cells=sampleLAD, n = "Cellular response to hypoxia in LAD", um=F, cex=0.5, logsc=T)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony, um=F, cex=0.5, logsc=T, cells=sampleLAD, "Hif1a", n = "Hif1a expression in LAD")
sampleCryo <- colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("MI_",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))]
plotfeatmap(sc1000_LAD_Cryo_Sham_harmony, Hypoxic_exp, cells=sampleCryo, n = "Cellular response to hypoxia in Cryoablation", um=F, cex=0.5, logsc=T)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony, um=F, cex=0.5, logsc=T, cells=sampleCryo, "Hif1a", n = "Hif1a expression in Cryoablation")
sampleSh <- colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)[grep("Sh_",colnames(sc1000_LAD_Cryo_Sham_harmony@ndata))]
plotfeatmap(sc1000_LAD_Cryo_Sham_harmony, Hypoxic_exp, cells=sampleSh, n = "Cellular response to hypoxia in Sham", um=F, cex=0.5, logsc=T)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony, um=F, cex=1, logsc=T, cells=sampleSh, "Hif1a", n = "Hif1a expression in Sham")


plotfeatmap(sc1000_LAD_Cryo_Sham_harmony, Hypoxic_exp, n = "Cellular response to hypoxia (LAD)", um=F, cex=0.2, logsc=T, cells = sampleLAD)
plotfeatmap(sc1000_LAD_Cryo_Sham_harmony, Hypoxic_exp, n = "Cellular response to hypoxia (Cryoablation)", um=F, cex=0.2, logsc=T, cells = sampleCryo)
plotfeatmap(sc1000_LAD_Cryo_Sham_harmony, Hypoxic_exp, n = "Cellular response to hypoxia (Sham)", um=F, cex=0.2, logsc=T, cells = sampleSh)


table(sc1000_LAD_Cryo_Sham_harmony@cpart)
violinMarkerPlot(Hypoxic_genes,sc1000_LAD_Cryo_Sham_harmony,set=c(1:31))



# plot normalized expressions across conditions
df1 <- data.frame(Hypoxic_exp[sampleLAD], "LAD")
colnames(df1) <- c("Hypoxic_gene_exp", "Injury_model")
df2 <- data.frame(Hypoxic_exp[sampleCryo], "Cryoablation")
colnames(df2) <- c("Hypoxic_gene_exp", "Injury_model")
df3 <- data.frame(Hypoxic_exp[sampleSh], "Sham")
colnames(df3) <- c("Hypoxic_gene_exp", "Injury_model")
df <- rbind(df1, df2, df3)

df$Injury_model <- as.factor(df$Injury_model)

# violin plot
p<-ggplot(df, aes(x=Injury_model, y=Hypoxic_gene_exp, fill=Injury_model)) +
  geom_violin(trim=T)+
  scale_x_discrete(limits=c("LAD", "Cryoablation", "Sham"))+
  geom_boxplot(width=0.4, fill="Grey")+
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 15),)+
  labs(title="Hypoxic genes expression",x="Injury model", y = "summed, normalized gene epxression")+
  ylim(0,0.02)
p + scale_fill_manual(values=c("#008BCC","#C9655E","#0A1722"))

# box plot
p<-ggplot(df, aes(x=Injury_model, y=Hypoxic_gene_exp, fill=Injury_model)) +
  scale_x_discrete(limits=c("LAD", "Cryoablation", "Sham"))+
  geom_boxplot(width=0.8)+
  theme(axis.text.x = element_text(vjust = 0.6,angle = 45, size = 15),)+
  labs(title="Hypoxic genes expression",x="Injury model", y = "summed, normalized gene epxression")+
  ylim(0,0.02)
p + scale_fill_manual(values=c("#008BCC","#C9655E","#0A1722"))





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








# plot gene expression heatmap 
require(pheatmap)
# create matrix (rows = injury model types, columns = gene names)
Hypoxic_genes_LADcells <- rowSums(sc1000_LAD_Cryo_Sham_harmony@ndata[Hypoxic_genes, sampleLAD]) / length(sampleLAD)
Hypoxic_genes_Cryocells <- rowSums(sc1000_LAD_Cryo_Sham_harmony@ndata[Hypoxic_genes, sampleCryo]) / length(sampleCryo)
Hypoxic_genes_Shamcells <- rowSums(sc1000_LAD_Cryo_Sham_harmony@ndata[Hypoxic_genes, sampleSh]) / length(sampleSh)
heatmap_matrix <- as.matrix(rbind(Hypoxic_genes_LADcells,Hypoxic_genes_Cryocells,Hypoxic_genes_Shamcells))
rownames(heatmap_matrix) <- c("LAD", "Cryoablation", "Sham")
pheatmap(heatmap_matrix, cluster_rows = FALSE, cluster_cols = T)
# plot in log scale
pheatmap(log(heatmap_matrix), cluster_rows = FALSE, cluster_cols = T)
# remove Ubb and Elob since it has much higher expression than other genes
pheatmap(log2(heatmap_matrix[,!colnames(heatmap_matrix)%in%c("Ubb", "Elob")]), cluster_rows = FALSE, cluster_cols = T)

# scale each gene's expression (score from 1 to -1), use "apply" and "scale" functions
# scale: scale each gene from -1 to 1
# apply: apply the scale function to all columns (2) in the matrix (ie. all genes)
heatmap_matrix_scale<-apply(heatmap_matrix,2,scale)
rownames(heatmap_matrix_scale) <- c("LAD", "Cryoablation", "Sham")
pheatmap(heatmap_matrix_scale, cluster_rows = FALSE, cluster_cols = T)
# scale by min-max (less pretty)
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
heatmap_matrix_minmax<-apply(heatmap_matrix,2,min_max_norm)
pheatmap(heatmap_matrix_minmax, cluster_rows = FALSE, cluster_cols = T)

#########################################################3


plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Pdgfra",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Pdgfra",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Ttn",logsc=T,fr=F, um=F, cex=0.2)
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_harmony,c(30),pvalue=0.01)
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_harmony,c(27),pvalue=0.01)

A <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)[sc1000_LAD_Cryo_Sham_harmony@cpart %in% c(2)]
B <- names(sc1000_LAD_Cryo_Sham_harmony@cpart)[sc1000_LAD_Cryo_Sham_harmony@cpart %in% c(7)]
x <- diffexpnb(getfdata(sc1000_LAD_Cryo_Sham_harmony,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Cl.2",Bname="Cl.7",show_names=TRUE,padj=TRUE)

# general tissue markers
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Ttn",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Pecam1",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Pdgfra",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Pdgfrb",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Acta2",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Ptprc",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Plp1",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Hbb-bt",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Postn",logsc=T,fr=F, um=F, cex=0.2)



# immune marker genes
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"S100a8",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Itgam",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Itgax",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Lyve1",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cx3cr1",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Retnla",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Itgae",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cd209a",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cd3e",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Il7r",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Il2ra",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Gata3",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Gzma",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cd19",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Mrc1",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Itgal",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Spp1",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Ccr2",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cd4",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cd8a",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Foxp3",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Trdc",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Hbb-bt",logsc=T,fr=F, um=F, cex=0.2)
plotexpmap(sc1000_LAD_Cryo_Sham_harmony,"Cpa3",logsc=T,fr=F, um=F, cex=0.2)



# cell type annotations

# CM
# 4, 27 (half)

# Fibroblasts
# 17,21,2,13,26,3,1,30

# Vascular cells
# 23,24,14,12,15

# neural cells
# 29

# immune cells general
# 22,6,20,8,18,5,26,7,10,19,28,16,30,25,11,9,31
# lymphocytes
# 8,13,15,16,20,24
# DC
# 17,5,21,19,11

annotated <- c(22,6,20,8,18,5,26,7,10,19,28,16,30,25,11,9,31,29,23,24,14,12,15,17,21,2,13,26,3,1,30,4,27)
annotated <- sort(annotated)
annotated 
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
# [24] 24 25 26 26 27 28 29 30 30 31
# 26 and 30 duplicated (26,30=FB/Immune)
table(sc1000_LAD_Cryo_Sham_harmony@cpart) # total 31 clusters



# immune cell clusters
cluster_ann_immuneALL <- c(cluster_ann_nonBlymphocytes, cluster_ann_DC, cluster_ann_immune)
cluster_ann_immuneALL <- cluster_ann_immuneALL[unique(names(cluster_ann_immuneALL))]
cluster_ann_immuneALL <- cluster_ann_immuneALL[names(sc1000_LAD_Cryo_Sham_harmony@cpart)]
unique(cluster_ann_immuneALL)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_immuneALL, um=T, cex=0.2, samples_col =  c(as.character(paletteer_c("grDevices::TealGrn", 30)), as.character(paletteer_c("grDevices::BluYl", 17)))
)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_immuneALL, um=T, cex=0.2, 
               subset  = c("T_Eff", "T_naive", "T_CD4_effector", "T_CD4_naive"))
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_immuneALL, um=T, cex=0.2, 
               subset  = c("MAIT", "MAIT_IL17"))

cluster_ann_immuneALL[cluster_ann_immuneALL%in%c("T_Eff", "T_naive")] <- "T_CD4_effector"
sort(unique(cluster_ann_immuneALL))



# rename CM subtypes
cluster_ann_CM1 <- cluster_ann_CM
names(cluster_ann_CM1) <- sub("Cr", "MI",names(cluster_ann_CM1))
names(cluster_ann_CM1) <- sub("CM_Ad_Sh", "CM1Ad_Sh",names(cluster_ann_CM1))
table(cluster_ann_CM1)
types <- str_sub(names(cluster_ann_CM1), 1,8)
unique(types)
cluster_ann_CM1[cluster_ann_CM1%in%"Angiogenic"] <- "CM_Angiogenic"
cluster_ann_CM1[cluster_ann_CM1%in%"Btk+"] <- "CM_Btk+"
cluster_ann_CM1[cluster_ann_CM1%in%"Dedifferentiating"] <- "CM_Dedifferentiating"
cluster_ann_CM1[cluster_ann_CM1%in%"Homeostatic"] <- "CM_Homeostatic"
cluster_ann_CM1[cluster_ann_CM1%in%"Hypertrophic"] <- "CM_Hypertrophic"
cluster_ann_CM1[cluster_ann_CM1%in%"Hypertrophic-pre"] <- "CM_Prehypertrophic"
cluster_ann_CM1[cluster_ann_CM1%in%"Slit2+"] <- "CM_Slit2+"


types <- str_sub(names(cluster_ann_CM1), 1,8)
unique(types)

# annotate cell types with all populations
cluster_ann <- sc1000_LAD_Cryo_Sham_harmony@cpart
length(cluster_ann)
types <- str_sub(names(cluster_ann), 1,8)
unique(types)


types <- str_sub(names(cluster_ann), 1,8)
unique(types)


cluster_ann <- c(cluster_ann_CM1, cluster_ann_Fibro, cluster_ann_neural, cluster_ann_Vascular, cluster_ann_immuneALL, cluster_ann)
cluster_ann <- cluster_ann[unique(names(cluster_ann))]
length(cluster_ann)
# in case cell number is not consistent, use this function to check
setdiff(names(cluster_ann),names(sc1000_LAD_Cryo_Sham_harmony@cpart))

types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 1,8)
unique(types)

cluster_ann <- cluster_ann[names(sc1000_LAD_Cryo_Sham_harmony@cpart)]
length(cluster_ann)
length(sc1000_LAD_Cryo_Sham_harmony@cpart)

cluster_ann_all <- cluster_ann
sort(unique(cluster_ann_all))
length(unique(names(cluster_ann)))
length(names(cluster_ann))
table(cluster_ann_all)

# manually annotated remaining cells:
dg <- clustdiffgenes(sc1000_LAD_Cryo_Sham_harmony,c(27),pvalue=0.01)

cluster_ann_all[cluster_ann_all%in%"4"] <- "CM_Homeostatic"

# annotate cluster 27 as dediff CM
#  cluster_ann_all[cl27_dediff_CM] <- "CM_Dedifferentiating"
#  cluster_ann_all[cl27_stressed] <- "Malat1_high"
# alternatively, keep cluster 27 not annotated
cluster_ann_all[cluster_ann_all%in%"27"] <- "Not_annotated"

# merge some annotations
cluster_ann_all[cluster_ann_all%in%c("DC_con1", "DC_con1_metabolic", "DC_con1_CD207")] <- "DC_con1"
cluster_ann_all[cluster_ann_all%in%c("DC_con2b", "DC_con2b_Ccl17", "DC_con2b_Retnla", "DC_IFN", "DC_mono", "DC_mono_Mgl2", "DC_neutrophil", "DC_Ccl4")] <- "DC_con2b/moDC"
cluster_ann_all[cluster_ann_all%in%c("FB_Macro_Saa3")] <- "FB_Saa3"
cluster_ann_all[cluster_ann_all%in%c("FB_Neutrophil", "FB_rMacro")] <- "Fibrocytes"
table(cluster_ann_all)

# create dataframe containing cell names, timepoints and annotations
cells_timepoints_annotations <- data.frame(cells_timepoints, cluster_ann_all)
head(cells_timepoints_annotations)
# export dataframe:
saveRDS(cells_timepoints_annotations, "mdata_cells_timepoints_annotations_detailed_1.rds")
# export VarID object:
saveRDS(sc1000_LAD_Cryo_Sham_harmony, "VarID_object_MI_allcells.rds")

# export subpopulations VarID object:
saveRDS(sc1000_LAD_Cryo_Sham_Immune1, "VarID_object_MI_Immune.rds")
saveRDS(sc1000_LAD_Cryo_Sham_nonBlymphocytes1, "VarID_object_MI_T_NK_ILC2.rds")
saveRDS(sc1000_LAD_Cryo_Sham_DCs, "VarID_object_MI_DC.rds")
saveRDS(sc1000_LAD_Cryo_Sham_Fibro, "VarID_object_MI_Fibroblasts.rds")
saveRDS(sc1000_LAD_Cryo_Sham_Vascular, "VarID_object_MI_Vascular.rds")
saveRDS(sc1000_LAD_Cryo_Sham_neural, "VarID_object_MI_SchwannCells.rds")
saveRDS(scCM_LAD_Cryo_sham1, "VarID_object_MI_CM.rds")





# plot the result:
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_all, um=F, cex=0.1,
               leg = T)
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_all, um=F, cex=0.1,
               subset = "CM_Ankrd1")

plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cells_timepoints_annotations$cluster_ann_all, um=F, cex=0.1,
               subset = "Neutrophil_2")

plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cells_timepoints_annotations$cluster_ann_all, um=F, cex=0.1,
               subset = "CM_Angiogenic")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_all, um=F, cex=0.1,
               subset = "CM_Dedifferentiating")
plotsymbolsmap(sc1000_LAD_Cryo_Sham_harmony, cluster_ann_all, um=F, cex=0.1,
               subset = "FB_transition_Il1rapl1")


# simplify cell types for NiCo analysis
table(cluster_ann_all)
table(cluster_ann_CM1)
table(cluster_ann_Fibro)
table(cluster_ann_neural)
table(cluster_ann_Vascular)
table(cluster_ann_immuneALL)


cluster_ann_simplify <- cluster_ann_all
cluster_ann_simplify[cluster_ann_simplify%in%c("CM_Homeostatic","CM_Angiogenic","CM_Slit2+","CM_Btk+")] <- "CM_General"
cluster_ann_simplify[cluster_ann_simplify%in%c("CM_Prehypertrophic", "CM_Hypertrophic")] <- "CM_BZ"
cluster_ann_simplify[cluster_ann_simplify%in%(names(table(cluster_ann_Fibro)))] <- "Fibroblasts"
cluster_ann_simplify[cluster_ann_simplify%in%(names(table(cluster_ann_neural)))] <- "SwC"
cluster_ann_simplify[cluster_ann_simplify%in%c("vEC_angio_IFN","vEC_Areg_Dkk2_Wnt","vEC_capillary1","vEC_capillary2")] <- "EC_capillary"
cluster_ann_simplify[cluster_ann_simplify%in%c("vPericyte_FB","vPericyte_INF","vPericyte_quiescent+")] <- "Pericytes"
cluster_ann_simplify[cluster_ann_simplify%in%c("vSMC_Ryr2","vSMC1","vSMC2")] <- "SMC"
cluster_ann_simplify[cluster_ann_simplify%in%c("Neutrophil_1","Neutrophil_2")] <- "Neutrophil"
cluster_ann_simplify[cluster_ann_simplify%in%c("mo_Macro_Il10_Fn1", "mo_Macro_Ly6C_Isg15_Cxcl3", "mo_Macro_Spp1")] <- "Macrophages_Ly6C"
cluster_ann_simplify[cluster_ann_simplify%in%c("Macro_Cxcl13","Macro_Retnla", "Macro_Timd4_Lyve1")] <- "Macrophages_Timd4"



cluster_ann_simplify[cluster_ann_simplify%in%c("NK_Gzma","NK_Klra5","NK_T")] <- "NK"
cluster_ann_simplify[cluster_ann_simplify%in%c("ILC2_IL5","ILC2")] <- "ILC2"
cluster_ann_simplify[cluster_ann_simplify%in%c("T_CD4_naive","T_CD4_effector")] <- "T_CD4"
cluster_ann_simplify[cluster_ann_simplify%in%c("T_CD8_naive","T_CD8_effector")] <- "T_CD8"
cluster_ann_simplify[cluster_ann_simplify%in%c("DC_Ccl4", "DC_con2_Ifitm1", "DC_con2b/moDC")] <- "DC_con2_mo"

table(cluster_ann_simplify)

# create dataframe containing cell names, timepoints and annotations (Simplified)
cells_timepoints_annotations_simplified <- data.frame(cells_timepoints, cluster_ann_simplify)
head(cells_timepoints_annotations_simplified)
# export dataframe:
saveRDS(cells_timepoints_annotations_simplified, "mdata_cells_timepoints_annotations_simplified_2.rds")




