require(CellChat)
# cell annotations
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


cell_id <- colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)
cell_id <- str_replace(cell_id, "MI", "Cr")
types <- str_sub(cell_id, 7,12)
unique(types)

# extract cells (and cluster annotations) of only a single time point
# examples of names: NM_Ad_MI_D56, CM1Ad_MI_D28, NM_Ad_Sh_D01
require(stringr)
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 7,12)
head(types)
unique(types)
cells_MI_D01 <- names(cluster_ann_all[types%in%c("MI_D01","LA_D01")])
length(cells_MI_D01)

# create a dataframe: rownames = cell names; column: cluster annotations
meta = data.frame(labels = cluster_ann_all, row.names = colnames(sc1000_LAD_Cryo_Sham_harmony@ndata)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels


# create Cellchat object of only Cryoablation+LAD D01 time point
# at the beginning (not after doing the step "subset data", the object name can only be "cellchat")
cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata[,cells_MI_D01], meta = meta[cells_MI_D01,,drop=F], group.by = "labels")
#> Create a CellChat object from a data matrix

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# optional:
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

cellchat_MI_D01 <- cellchat

# repeat CellChat with Cryoablation+LAD D03
cells_MI_D03 <- names(cluster_ann_all[types%in%c("MI_D03", "LA_D03")])
length(cells_MI_D03)

cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata[,cells_MI_D03], meta = meta[cells_MI_D03,,drop=F], group.by = "labels")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)

cellchat_MI_D03 <- cellchat

# repeat CellChat with Cryoablation+LAD D07
cells_MI_D07 <- names(cluster_ann_all[types%in%c("MI_D07", "LA_D07")])
length(cells_MI_D07)

cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata[,cells_MI_D07], meta = meta[cells_MI_D07,,drop=F], group.by = "labels")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)

cellchat_MI_D07 <- cellchat

# repeat CellChat with Cryoablation+LAD D28
cells_MI_D28 <- names(cluster_ann_all[types%in%c("MI_D28", "LA_D28")])
length(cells_MI_D28)

cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata[,cells_MI_D28], meta = meta[cells_MI_D28,,drop=F], group.by = "labels")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)

cellchat_MI_D28 <- cellchat

# repeat CellChat with Cryoablation+LAD D56
cells_MI_D56 <- names(cluster_ann_all[types%in%c("MI_D56", "LA_D56")])
length(cells_MI_D56)

cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata[,cells_MI_D56], meta = meta[cells_MI_D56,,drop=F], group.by = "labels")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)

cellchat_MI_D56 <- cellchat

# repeat CellChat with Sham
types <- str_sub(colnames(sc1000_LAD_Cryo_Sham_harmony@ndata), 7,8)
head(types)
unique(types)
cells_Sh <- names(cluster_ann_all[types%in%"Sh"])
length(cells_Sh)

cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata[,cells_Sh], meta = meta[cells_Sh,,drop=F], group.by = "labels")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)

cellchat_Sh <- cellchat


# merging all Cellchat objects into one
cellchat <- createCellChat(object = sc1000_LAD_Cryo_Sham_harmony@ndata, meta = meta, group.by = "labels")

# Define the cell labels to lift up
group.new = levels(cellchat@idents)
cellchat_MI_D01 <- liftCellChat(cellchat_MI_D01, group.new)
cellchat_MI_D03 <- liftCellChat(cellchat_MI_D03, group.new)
cellchat_MI_D07 <- liftCellChat(cellchat_MI_D07, group.new)
cellchat_MI_D28 <- liftCellChat(cellchat_MI_D28, group.new)
cellchat_MI_D56 <- liftCellChat(cellchat_MI_D56, group.new)
cellchat_Sh <- liftCellChat(cellchat_Sh, group.new)

cellchat_MI_D01 <- netAnalysis_computeCentrality(cellchat_MI_D01)
cellchat_MI_D03 <- netAnalysis_computeCentrality(cellchat_MI_D03)
cellchat_MI_D07<- netAnalysis_computeCentrality(cellchat_MI_D07)
cellchat_MI_D28 <- netAnalysis_computeCentrality(cellchat_MI_D28)
cellchat_MI_D56 <- netAnalysis_computeCentrality(cellchat_MI_D56)
cellchat_Sh <- netAnalysis_computeCentrality(cellchat_Sh)


object.list <- list(MI_D01 = cellchat_MI_D01, MI_D03 = cellchat_MI_D03, MI_D07 = cellchat_MI_D07, MI_D28 = cellchat_MI_D28, MI_D56 = cellchat_MI_D56, Sham = cellchat_Sh)
# always run "object.list" AFTER computing centrality of individual cellchat objects!


cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)



gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1:6))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1:6), measure = "weight")
gg1 + gg2




# Plot number and strength of interactions
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# show all cell types
sort(unique(cellchat@meta$labels))
names(object.list)

# Neural outward signals

netVisual_circle(object.list[[1]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[1]),top = 0.1, sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[2]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[2]),top = 0.1, sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[3]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[3]),top = 0.1, sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[4]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[4]),top = 0.1, sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[5]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[5]),top = 0.1, sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[6]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[6]),top = 0.1, sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)

# Neural-to-Lymphocyte signals

netVisual_circle(object.list[[1]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[1]), idents.use = c("T_CD4_effector","T_CD8_naive", "T_CD4_naive", "T_CD8_effector", "T_IFNg_naive", "T_gd", "T_Isg15", "Treg", "NK_T", "NK_Gzma", "NK_Klra5", "ILC2", "ILC2_IL5", "MAIT", "MAIT_IL17"), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[2]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[2]), idents.use = c("T_CD4_effector","T_CD8_naive", "T_CD4_naive", "T_CD8_effector", "T_IFNg_naive", "T_gd", "T_Isg15", "Treg", "NK_T", "NK_Gzma", "NK_Klra5", "ILC2", "ILC2_IL5", "MAIT", "MAIT_IL17"), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[3]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[3]), idents.use = c("T_CD4_effector","T_CD8_naive", "T_CD4_naive", "T_CD8_effector", "T_IFNg_naive", "T_gd", "T_Isg15", "Treg", "NK_T", "NK_Gzma", "NK_Klra5", "ILC2", "ILC2_IL5", "MAIT", "MAIT_IL17"), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[4]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[4]), idents.use = c("T_CD4_effector","T_CD8_naive", "T_CD4_naive", "T_CD8_effector", "T_IFNg_naive", "T_gd", "T_Isg15", "Treg", "NK_T", "NK_Gzma", "NK_Klra5", "ILC2", "ILC2_IL5", "MAIT", "MAIT_IL17"), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[5]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[5]), idents.use = c("T_CD4_effector","T_CD8_naive", "T_CD4_naive", "T_CD8_effector", "T_IFNg_naive", "T_gd", "T_Isg15", "Treg", "NK_T", "NK_Gzma", "NK_Klra5", "ILC2", "ILC2_IL5", "MAIT", "MAIT_IL17"), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)
netVisual_circle(object.list[[6]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[6]), idents.use = c("T_CD4_effector","T_CD8_naive", "T_CD4_naive", "T_CD8_effector", "T_IFNg_naive", "T_gd", "T_Isg15", "Treg", "NK_T", "NK_Gzma", "NK_Klra5", "ILC2", "ILC2_IL5", "MAIT", "MAIT_IL17"), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"), remove.isolate = T)

# Neural-to-CM signals

netVisual_circle(object.list[[1]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[1]), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"),targets.use = c("CM_Angiogenic","CM_Ankrd1","CM_Dedifferentiating","CM_Gck","CM_Hypertrophic","CM_IFN","CM_Metabolic","CM_Myh7","CM_Normal","CM_Ppargc1a","CM_Rasef","CM_Slit2"), remove.isolate = T)
netVisual_circle(object.list[[2]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[2]), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"),targets.use = c("CM_Angiogenic","CM_Ankrd1","CM_Dedifferentiating","CM_Gck","CM_Hypertrophic","CM_IFN","CM_Metabolic","CM_Myh7","CM_Normal","CM_Ppargc1a","CM_Rasef","CM_Slit2"), remove.isolate = T)
netVisual_circle(object.list[[3]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[3]), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"),targets.use = c("CM_Angiogenic","CM_Ankrd1","CM_Dedifferentiating","CM_Gck","CM_Hypertrophic","CM_IFN","CM_Metabolic","CM_Myh7","CM_Normal","CM_Ppargc1a","CM_Rasef","CM_Slit2"), remove.isolate = T)
netVisual_circle(object.list[[4]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[4]), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"),targets.use = c("CM_Angiogenic","CM_Ankrd1","CM_Dedifferentiating","CM_Gck","CM_Hypertrophic","CM_IFN","CM_Metabolic","CM_Myh7","CM_Normal","CM_Ppargc1a","CM_Rasef","CM_Slit2"), remove.isolate = T)
netVisual_circle(object.list[[5]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[5]), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"),targets.use = c("CM_Angiogenic","CM_Ankrd1","CM_Dedifferentiating","CM_Gck","CM_Hypertrophic","CM_IFN","CM_Metabolic","CM_Myh7","CM_Normal","CM_Ppargc1a","CM_Rasef","CM_Slit2"), remove.isolate = T)
netVisual_circle(object.list[[6]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12,vertex.label.cex = 0.5, title.name = paste0("Number of interactions - ", names(object.list)[6]), sources.use =c("Schwann_Galectin", "Schwann_IFN", "Schwann_metabolic", "Schwann_quiescent"),targets.use = c("CM_Angiogenic","CM_Ankrd1","CM_Dedifferentiating","CM_Gck","CM_Hypertrophic","CM_IFN","CM_Metabolic","CM_Myh7","CM_Normal","CM_Ppargc1a","CM_Rasef","CM_Slit2"), remove.isolate = T)



# Compare outgoing (or incoming) signaling associated with each cell population

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
# pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# outgoing signals
ht0 = netAnalysis_signalingRole_heatmap(object.list[[i+5]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+5], width = 8, height = 18, font.size = 5)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 18, font.size = 5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1],width = 8, height = 18, font.size = 5)
draw(ht0 + ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 8, height = 18, font.size = 5)
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 8, height = 18, font.size = 5)
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+4], width = 8, height = 18, font.size = 5)
draw(ht3 + ht4 + ht5, ht_gap = unit(0.5, "cm"))

# incoming signals
ht0 = netAnalysis_signalingRole_heatmap(object.list[[i+5]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+5], width = 8, height = 18, font.size = 5, color.heatmap = "GnBu")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 18, font.size = 5, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 18, font.size = 5, color.heatmap = "GnBu")
draw(ht0 + ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 8, height = 18, font.size = 5, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 8, height = 18, font.size = 5, color.heatmap = "GnBu")
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+4], width = 8, height = 18, font.size = 5, color.heatmap = "GnBu")
draw(ht3 + ht4 + ht5, ht_gap = unit(0.5, "cm"))




# list all cell types
sort(unique(cellchat@meta$labels))

# Identify dysfunctional signaling by comparing the communication probabities
sort(unique(cellchat@meta$labels))
names(object.list)

# show all pathways
sort(unique(c(cellchat@netP$MI_D01$pathways, cellchat@netP$MI_D03$pathways, cellchat@netP$MI_D07$pathways, cellchat@netP$MI_D28$pathways, cellchat@netP$MI_D56$pathways, cellchat@netP$Sham$pathways)))


Neural <- sort(unique(cellchat@meta$labels))[66:69]
FB <- sort(unique(cellchat@meta$labels))[29:47]
CM <- sort(unique(cellchat@meta$labels))[2:13]
DC <- sort(unique(cellchat@meta$labels))[14:27]
Macrophages <- sort(unique(cellchat@meta$labels))[c(50:54,58:60)]
Neutrophils <- sort(unique(cellchat@meta$labels))[61:62]
Mast_cells <- sort(unique(cellchat@meta$labels))[57]
Myeloids <- c(Macrophages, Neutrophils, DC, Mast_cells)
Lymphocytes <- sort(unique(cellchat@meta$labels))[c(1, 48,49, 55,56, 63:65, 70:78)]
Vascular <- sort(unique(cellchat@meta$labels))[79:94]





# plot cell type interaction dotplots

# Macrophages/Neutrophils --> hypertrophic CM (Day 1, 3 sham)
netVisual_bubble(cellchat, sources.use = c(Myeloids), targets.use = "CM_Hypertrophic",  comparison = c(1,2,6), angle.x = 45,font.size = 8)
netVisual_bubble(cellchat, sources.use = c(Macrophages, Mast_cells, Neutrophils), targets.use = "CM_Hypertrophic",  comparison = c(1,2,6), angle.x = 45,font.size = 12)
netVisual_bubble(cellchat, sources.use = c(FB), targets.use = "CM_Hypertrophic",  comparison = c(1,2,6), angle.x = 45,font.size = 4)


# FB --> Macrophages
netVisual_bubble(cellchat, sources.use = c("FBmyo","FB_transition_Cd9","FB_quiescent"), targets.use = c("mo_Macro_Ly6C_Isg15_Cxcl3","Macro_Trem2_Gpnmb", "Macro_MHCII_Cx3cr1", "Macro_Timd4_Lyve1"),  comparison = c(3,4), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, sources.use = c("FBmyo","FB_transition_Cd9","FB_quiescent"), targets.use = c("mo_Macro_Ly6C_Isg15_Cxcl3","Macro_Trem2_Gpnmb", "Macro_MHCII_Cx3cr1", "Macro_Timd4_Lyve1"),  comparison = c(2,3,4), angle.x = 45,font.size = 6, signaling = c("FGF", "SEMA3", "SEMA4", "BMP", "ANGPT", "FN1"))
netVisual_bubble(cellchat, sources.use = c("FBmyo","FB_transition_Cd9","FB_quiescent"), targets.use = c("mo_Macro_Ly6C_Isg15_Cxcl3","Macro_Trem2_Gpnmb", "Macro_MHCII_Cx3cr1", "Macro_Timd4_Lyve1"),  comparison = c(1,2), angle.x = 45,font.size = 6, signaling = c("FGF", "SEMA3", "SEMA4", "BMP", "ANGPT", "FN1"))
# Macrophages --> FB
netVisual_bubble(cellchat, targets.use = c(FB), sources.use = c(Macrophages),  comparison = c(1,2), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, targets.use = c("FBmyo","FB_transition_Cd9","FB_quiescent"), sources.use = c("mo_Macro_Ly6C_Isg15_Cxcl3","Macro_Trem2_Gpnmb", "Macro_MHCII_Cx3cr1", "Macro_Timd4_Lyve1", Neutrophils),  comparison = c(1,2,6), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, targets.use = c("FBmyo","FB_transition_Cd9","FB_quiescent"), sources.use = c("mo_Macro_Ly6C_Isg15_Cxcl3","Macro_Trem2_Gpnmb", "Macro_MHCII_Cx3cr1", "Macro_Timd4_Lyve1"),  comparison = c(2,3,4), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, targets.use = c("FBmyo","FB_transition_Cd9","FB_quiescent"), sources.use = c("mo_Macro_Ly6C_Isg15_Cxcl3","Macro_Trem2_Gpnmb", "Macro_MHCII_Cx3cr1", "Macro_Timd4_Lyve1"),  comparison = c(2,3,4), angle.x = 45,font.size = 6, signaling = c("GAS", "PROS"))


# MAIT --> others
netVisual_bubble(cellchat, sources.use = c("MAIT", "MAIT_IL17"),  comparison = c(1,2), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat,  targets.use = DC, sources.use = c("MAIT", "MAIT_IL17"),  comparison = c(1:6), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat,  targets.use = c(DC, Lymphocytes), sources.use = c("MAIT_IL17"),  comparison = c(3,6), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat,  targets.use = c(CM), sources.use = c("MAIT", "MAIT_IL17"),  comparison = c(3,6), angle.x = 45,font.size = 6)



# SwC --> lymphocytes
netVisual_bubble(cellchat, sources.use = Neural, targets.use = Lymphocytes,  comparison = c(1,6), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("ILC2", "ILC2_IL5"),  comparison = c(1:6), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("MAIT", "MAIT_IL17"),  comparison = c(1:6), angle.x = 45,font.size = 6)

# SwC --> CM
netVisual_bubble(cellchat, sources.use = Neural, targets.use = CM,  comparison = c(1,2,6), angle.x = 45,font.size = 4)
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("CM_Dedifferentiating", "CM_Hypertrophic", "CM_Normal"),  comparison = c(1,2,6), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("CM_Dedifferentiating","CM_Hypertrophic", "CM_Normal"),  comparison = c(3,4,6), angle.x = 45,font.size = 6)
#plot specific pathways only
# check which pathway of specified ligand/receptor belongs to (use capital letters)
head(cellchat@DB$interaction[grep("BTC",cellchat@DB$interaction$interaction_name),"pathway_name"])
head(cellchat@DB$interaction[grep("WNT6",cellchat@DB$interaction$interaction_name),"pathway_name"])
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("CM_Hypertrophic", "CM_Normal"),  comparison = c(1:6), angle.x = 45,font.size = 6, signaling = c("EGF", "WNT"))

# SwC --> Vascular cells
netVisual_bubble(cellchat, sources.use = Neural, targets.use = Vascular,  comparison = c(1,6), angle.x = 45,font.size = 4)
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("vEC_angio_IFN","vEC_Areg_Dkk2_Wnt", "vEC_Arterial", "vEC_capillary1"),  comparison = c(1:6), angle.x = 45,font.size = 4, signaling = c("EGF", "WNT"))
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("vEC_capillary2","vEC_Endocardial", "vEC_Lymphatic"),  comparison = c(1:6), angle.x = 45,font.size = 4, signaling = c("EGF", "WNT"))
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("vEC_metabolic","vEpicardial_derived", "vEC_Immune"),  comparison = c(1:6), angle.x = 45,font.size = 4, signaling = c("EGF", "WNT"))
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("vPericyte_FB","vPericyte_INF", "vPericyte_quiescent"),  comparison = c(1:6), angle.x = 45,font.size = 4, signaling = c("EGF", "WNT"))
netVisual_bubble(cellchat, sources.use = Neural, targets.use = c("vSMC_Ryr2","vSMC1","vSMC2"),  comparison = c(1:6), angle.x = 45,font.size = 4, signaling = c("EGF", "WNT"))
# shortlisted plots
netVisual_bubble(cellchat, sources.use = c("Schwann_Galectin"), targets.use = c("vPericyte_FB","vPericyte_INF", "vPericyte_quiescent"),  comparison = c(1:6), angle.x = 45,font.size = 10, signaling = c("EGF"))
netVisual_bubble(cellchat, sources.use = c("Schwann_metabolic", "Schwann_quiescent"), targets.use = c("vEC_Areg_Dkk2_Wnt", "vEC_Arterial", "vEC_capillary1","vEC_Lymphatic","vSMC2"),  comparison = c(1:6), angle.x = 45,font.size = 6, signaling = c("WNT"))



netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("Fibroblast"),  comparison = c(2:4), angle.x = 45,font.size = 6, signaling = c("TENASCIN","ncWNT","SELL","GAS","PROS","NEGR","MK"))
netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("CM_general","CM_hypertrophic"),  comparison = c(3), angle.x = 45,font.size = 6)

netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("monocyte", "Macrophage_Cx3cr1", "Macrophage_Lyve1", "Macrophage_Retnla", "Macrophage_DC_cycling"),  comparison = c(3), angle.x = 45,font.size = 6, signaling = c("CD45", "TENASCIN", "CSF", "CHEMERIN", "NOTCH"))
netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("Macrophage_Cx3cr1", "Macrophage_Lyve1"),  comparison = c(4), angle.x = 45,font.size = 8, color.text.use = F, signaling = c("CD45", "TENASCIN", "CSF", "CHEMERIN", "NOTCH", "COMPLEMENT", "GDF",  "APRIL", "IL4"))
netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("monocyte", "Macrophage_Lyve1"),  comparison = c(3), angle.x = 45,font.size = 8, color.text.use = F, signaling = c("CD45", "TENASCIN", "CSF", "CHEMERIN", "NOTCH", "COMPLEMENT", "GDF",  "APRIL", "IL4"))
netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("monocyte", "Macrophage_Lyve1"),  comparison = c(4), angle.x = 45,font.size = 8, color.text.use = F, signaling = c("CD45", "TENASCIN", "CSF", "CHEMERIN", "NOTCH", "COMPLEMENT", "GDF",  "APRIL", "IL4"))
netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("CM_dediff", "CM_hypertrophic"),  comparison = c(1,3), angle.x = 45,font.size = 7, color.text.use = T, signaling = c("FN1","JAM","NCAM","VTN","TWEAK","MK","VCAM","PERIOSTIN","PTN"))

netVisual_bubble(cellchat, sources.use = c("Fibroblast", "Fibroblast_Saa","Neural"), targets.use = c("monocyte", "Macrophage_DC_cycling", "Macrophage_Lyve1", "Macrophage_Cx3cr1", "Macrophage_Retnla"),  comparison = c(3:4), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, sources.use = c("Fibroblast", "Fibroblast_Saa","Neural"), targets.use = c("monocyte", "Macrophage_DC_cycling", "Macrophage_Lyve1", "Macrophage_Cx3cr1", "Macrophage_Retnla"),  comparison = c(3:4), angle.x = 45,font.size = 6, signaling = c("APP", "CD200", "SEMA3", "GAS", "FGF"))

netVisual_bubble(cellchat, sources.use = c("T_gd"), targets.use = c(1:33), comparison = c(3, 5), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, targets.use = c("T_gd"), sources.use = c(1:33), comparison = c(3, 5), angle.x = 45,font.size = 6)

netVisual_bubble(cellchat, sources.use = c("ILC2"),  comparison = c(3), angle.x = 45,font.size = 6)
netVisual_bubble(cellchat, sources.use = c("ILC2"), targets.use = Macrophages, comparison = c(2,3,4), angle.x = 45,font.size = 6)

netVisual_bubble(cellchat, sources.use = c("ILC2", "MAIT", "MAIT_IL17", "T_gd", "Treg"), targets.use = "Macro_MHCII_Cx3cr1",  comparison = c(2,3,4), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, sources.use = c("ILC2", "MAIT", "MAIT_IL17", "T_gd", "Treg"), targets.use = "Macro_Timd4_Lyve1",  comparison = c(2,3,4), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, targets.use = c("DC_con1","DC_CCR7","ILC2", "MAIT", "MAIT_IL17", "Treg"), sources.use = "vEC_Lymphatic",  comparison = c(2,3,4), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, targets.use = c("ILC2", "MAIT", "MAIT_IL17", "T_gd", "Treg"), sources.use = "DC_CCR7",  comparison = c(2,3,4), angle.x = 45,font.size = 10)


netVisual_bubble(cellchat, targets.use = c("CM_Dedifferentiating", "CM_Hypertrophic"), sources.use = c("ILC2", "ILC2_IL5", "T_gd"),  comparison = c(2,3,4), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, targets.use = c("CM_Dedifferentiating"), sources.use = c("Schwann_Galectin", "Schwann_IFN", "Schwann_quiescent"),  comparison = c(2,3,4), angle.x = 45,font.size = 10)


netVisual_bubble(cellchat, sources.use = c("monocyte","EC_micro", "NK_cell"), targets.use = c("CM_general","CM_hypertrophic"),  comparison = c(1:6), angle.x = 45,font.size = 6, signaling = "GALECTIN" )

netVisual_bubble(cellchat, sources.use = c("Fibroblast", "Fibroblast_Saa","Neural"), targets.use = c("monocyte", "Macrophage_DC_cycling", "Macrophage_Lyve1", "Macrophage_Cx3cr1", "Macrophage_Retnla"),  comparison = c(3:4), angle.x = 45,font.size = 6, signaling = c("APP", "CD200", "SEMA3", "GAS", "FGF"))

netVisual_bubble(cellchat, sources.use = c("ILC2"), targets.use = c("DC_Ccr7", "cDC_moDC", "cDC1"),  comparison = c(1:3), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, sources.use = c("ILC2"), targets.use = c("CM_general","CM_hypertrophic"),  comparison = c(1:6), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, sources.use = c("ILC2"), targets.use = c("Fibroblast_Saa","Fibroblast"),  comparison = c(1:6), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, targets.use = c("ILC2", "T_gd"), sources.use = c("Fibroblast"),  comparison = c(1:6), angle.x = 45,font.size = 10)
netVisual_bubble(cellchat, sources.use = c("ILC2"), targets.use = c("monocyte","Macrophage_Lyve1", "Macrophage_Cx3cr1"),  comparison = c(1:6), angle.x = 45,font.size = 10)

netVisual_bubble(cellchat, sources.use = c("cDC_moDC","cDC1"), targets.use = c("T_gd","Treg", "ILC2", "T_helper", "T_killer", "T_naive_memory", "T_cell_cycling"),  comparison = c(3), angle.x = 45,font.size = 10)

netVisual_bubble(cellchat, sources.use = c(1:33), targets.use = c("Fibroblast"), comparison = c(2,4), angle.x = 45, font.size = 6, signaling = c("EPHB", "CHEMERIN", "OSM", "GDF", "CD46", "EDN", "ACTIVIN"))
netVisual_bubble(cellchat, sources.use = c("Macrophage_Lyve1","SMC", "Pericyte", "CM_hypertrophic", "CM_Ppargc1a", "CM_general", "CM_dediff", "EC_macro"), targets.use = c("Fibroblast"), comparison = c(1:6), angle.x = 45, font.size = 6, signaling = c("EPHB", "CHEMERIN", "OSM", "GDF", "CD46", "EDN", "ACTIVIN", "TENASCIN","ncWNT","GAS","PROS"))



netVisual_bubble(cellchat, sources.use = c("FB_quiescent"), targets.use = c("Macro_MHCII_Cx3cr1"),  comparison = c(3,4), angle.x = 45,font.size = 6)

netVisual_bubble(cellchat, sources.use = c("FBmyo"), targets.use = c("Macro_MHCII_Cx3cr1"),  comparison = c(3,4), angle.x = 45,font.size = 6)

netVisual_bubble(cellchat, sources.use = c("FB_quiescent"), targets.use = c("Macro_MHCII_Cx3cr1"),  comparison = c(3,4), angle.x = 45,font.size = 6)

netVisual_bubble(cellchat, sources.use = FB, targets.use = Macrophages,  comparison = c(3,4), angle.x = 45,font.size = 6, signaling = c("GAS","PROS"))


