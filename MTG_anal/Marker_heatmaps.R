if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

library(scrattch.mapping)

#X <-matrix(rpois(100,10),ncol=10)
#colnames(X) = c("A","B","C","D","E","F","G","H","I","J")
load(file.path(mappingFolder,"NHP_BG_204_329_AIT115ann_map2.Rdata"))

anno_conf = anno_mapped_sub[anno_mapped_sub$level3.subclass_Corr =='D1-Matrix' &
                              anno_mapped_sub$level3.subclass_Tree == 'D2-Matrix',]
anno_conf$source = "D1/D2 conf samples"
#rownames(anno_conf)
load(paste0(data_dir, "/20230518_RSC-204-329_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20230518_RSC-204-329_macaque_patchseq_star2.7_samp.dat.Rdata"))

query.metadata <- samp.dat
counts      <- cpmR   # Genes are rows, samples are columns

query.counts   <- counts
query.data   <- logCPM(query.counts)

# Put annotations and counts in the same order
query.metadata <- query.metadata[match(colnames(query.data),query.metadata$exp_component_name),] 
rownames(query.metadata) <- query.metadata$exp_component_name 
query.data<-t(query.data)
gene_list = c("DRD1", "DRD2", "TAC1", "PENK", "TAC3", "RXFP1", "CPNE4", "STXBP6", "KCNIP1")
data_conf <- query.data[anno_conf$exp_component_name,gene_list]
data_conf <- as.data.frame(data_conf)
data_conf$source = "D1/D2 conf samples"

AIT.anndata <- loadTaxonomy(refFolder)
D1keep = which(AIT.anndata$obs$level3.subclass_label == "D1-Matrix")
D1keep = sample(D1keep, 35)
anndata_D1 <- AIT.anndata$X[D1keep,gene_list]
anndata_D1 <- as.data.frame(anndata_D1)
anndata_D1$source = "D1-Matrix taxonomy"

D2keep = which(AIT.anndata$obs$level3.subclass_label == "D2-Matrix")
D2keep = sample(D2keep, 35)
anndata_D2 <- AIT.anndata$X[D2keep,gene_list]
anndata_D2 <- as.data.frame(anndata_D2)
anndata_D2$source = "D2-Matrix taxonomy"

D1D2keep = which(AIT.anndata$obs$level3.subclass_label == "D1D2-Hybrid")
D1D2keep = sample(D1D2keep, 35)
anndata_D1D2 <- AIT.anndata$X[D1D2keep,gene_list]
anndata_D1D2 <- as.data.frame(anndata_D1D2)
anndata_D1D2$source = "D1D2-Hybrid taxonomy"

all <- rbind(data_conf, anndata_D1, anndata_D2, anndata_D1D2)


#colmap <- colorRamp2(c(min(X), 0, max(X)), c("blue", "white", "red"))
colmap <- colorRamp2(c(min(all[,1:9]), max(all[,1:9])), c("white", "red"))
#ha <- HeatmapAnnotation()

#h1 <- Heatmap(X, name = "test", cluster_columns = T, column_labels = colnames(X), col = colmap)
#cluster_rows
# or use default colormap

jpeg(file.path(mappingFolder,'204_329_D1corr_D2tree_heatmap.jpg'), quality = 100, width = 750, height = 750)

Heatmap(all[,1:9], cluster_rows = T, cluster_columns = F, cluster_row_slices = T, column_labels = colnames(all)[1:9], show_row_names = FALSE, col = colmap, row_split = all[,10])

dev.off()


# scProjection conflicts
gene_list = c("PPP1R1B", "BCL11B", "PDE1B", "DRD2", "PENK", "DRD1", "TAC1", "RXFP1", "CPNE4", "STXBP6", "KCNIP1")
sc_conf_cells = c("SM-J2NRM_S201_E1-50",
              "SM-J3A23_S154_E1-50",
              "SM-IRUWE_S706_E1-50",
              "SM-IRUWE_S707_E1-50",
              "SM-IRUWE_S708_E1-50",
              "SM-IRUWE_S767_E1-50",
              "AB-S40302_S472_E1-50",
              "AB-S40302_S477_E1-50",
              "AB-S40302_S478_E1-50",
              "AB-S40302_S604_E1-50",
              "AB-S40302_S605_E1-50",
              "AB-S40303_S724_E1-50",
              "AB-S40303_S728_E1-50",
              "AB-S40303_S730_E1-50",
              "AB-S40304_S341_E1-50",
              "AB-S40305_S118_E1-50",
              "AB-S40305_S132_E1-50",
              "AB-S40305_S135_E1-50",
              "AB-S40305_S708_E1-50",
              "AB-S40306_S597_E1-50",
              "AB-S40306_S668_E1-50",
              "AB-S40306_S669_E1-50",
              "AB-S40308_S256_E1-50",
              "AB-S40308_S466_E1-50",
              "AB-S40310_S721_E1-50",
              "AB-S40311_S343_E1-50",
              "AB-S40312_S445_E1-50",
              "AB-S40312_S446_E1-50",
              "AB-S40312_S449_E1-50",
              "AB-S40318_S433_E1-50") 

morpho_conf_cells = c("AB-S40311_S595_E1-50")

data_conf <- query.data[sc_conf_cells,gene_list]
data_conf <- as.data.frame(data_conf)
data_conf$source = "scProjection/Corr conflicts"

data_morpho_conf <- query.data[morpho_conf_cells,gene_list]
data_morpho_conf <- data.frame(t(as.data.frame(data_morpho_conf)))  # To get correct orientation for vector
data_morpho_conf$source = "Morpho/Scrattch conflict"

D1keep = which(AIT.anndata$obs$level3.subclass_label == "D1-Matrix")
D1keep = sample(D1keep, 35)
anndata_D1 <- AIT.anndata$X[D1keep,gene_list]
anndata_D1 <- as.data.frame(anndata_D1)
anndata_D1$source = "D1-Matrix taxonomy"

D2keep = which(AIT.anndata$obs$level3.subclass_label == "D2-Matrix")
D2keep = sample(D2keep, 35)
anndata_D2 <- AIT.anndata$X[D2keep,gene_list]
anndata_D2 <- as.data.frame(anndata_D2)
anndata_D2$source = "D2-Matrix taxonomy"

#D1D2keep = which(AIT.anndata$obs$level3.subclass_label == "D1D2-Hybrid")
#D1D2keep = sample(D1D2keep, 35)
#anndata_D1D2 <- AIT.anndata$X[D1D2keep,gene_list]
#anndata_D1D2 <- as.data.frame(anndata_D1D2)
#anndata_D1D2$source = "D1D2-Hybrid taxonomy"

INkeep = which(AIT.anndata$obs$level1.class_label == "IN")
INkeep = sample(INkeep, 50)
anndata_IN <- AIT.anndata$X[INkeep,gene_list]
anndata_IN <- as.data.frame(anndata_IN)
anndata_IN$source = "Mixed IN"

all <- rbind(data_conf, data_morpho_conf, anndata_D1, anndata_D2, anndata_IN)

colmap <- colorRamp2(c(min(all[,1:11]), max(all[,1:11])), c("white", "red"))

jpeg(file.path(mappingFolder,'204_329_MSN_IN_conflicts.jpg'), quality = 100, width = 750, height = 750)

par(mar = c(5.1, 4.1, 7.1, 2.1))

ht = Heatmap(all[,1:11], cluster_rows = T, cluster_columns = F, cluster_row_slices = T, 
             column_labels = colnames(all)[1:11], show_row_names = FALSE, col = colmap, 
             row_split = all[,12], row_gap=unit(.05, "npc"))
draw(ht, padding = unit(c(2, 2, 25, 2), "mm"))

dev.off()


# With Pvalb-col19A1 genes and samples to compare for the morpho conflict sample
gene_list = c("PPP1R1B", "BCL11B", "PDE1B", "DRD2", "PENK", "DRD1", "TAC1", "RXFP1", 
"CPNE4", "STXBP6", "KCNIP1","PVALB", "COL19A1", "ST18")

data_morpho_conf <- query.data[morpho_conf_cells,gene_list]
data_morpho_conf <- data.frame(t(as.data.frame(data_morpho_conf)))  # To get correct orientation for vector
data_morpho_conf$source = "Cell #1199243739"

D1keep = which(AIT.anndata$obs$level3.subclass_label == "D1-Matrix")
D1keep = sample(D1keep, 35)
anndata_D1 <- AIT.anndata$X[D1keep,gene_list]
anndata_D1 <- as.data.frame(anndata_D1)
anndata_D1$source = "D1-Matrix taxonomy"

D2keep = which(AIT.anndata$obs$level3.subclass_label == "D2-Matrix")
D2keep = sample(D2keep, 35)
anndata_D2 <- AIT.anndata$X[D2keep,gene_list]
anndata_D2 <- as.data.frame(anndata_D2)
anndata_D2$source = "D2-Matrix taxonomy"

INkeep = which(AIT.anndata$obs$level1.class_label == "IN")
INkeep = sample(INkeep, 50)
anndata_IN <- AIT.anndata$X[INkeep,gene_list]
anndata_IN <- as.data.frame(anndata_IN)
anndata_IN$source = "Mixed IN"

PVALBkeep = which(AIT.anndata$obs$level3.subclass_label == "PVALB-COL19A1-ST18")
PVALBkeep = sample(PVALBkeep, 50)
anndata_PVALB <- AIT.anndata$X[PVALBkeep,gene_list]
anndata_PVALB <- as.data.frame(anndata_PVALB)
anndata_PVALB$source = "PVALB-COL19A1-ST18"

all <- rbind(data_morpho_conf, anndata_D1, anndata_D2, anndata_IN, anndata_PVALB)

colmap <- colorRamp2(c(min(all[,1:14]), max(all[,1:14])), c("white", "red"))

jpeg(file.path(mappingFolder,'204_329_morpho_conflict.jpg'), quality = 100, width = 750, height = 750)

par(mar = c(5.1, 4.1, 7.1, 2.1))

ht = Heatmap(all[,1:14], cluster_rows = T, cluster_columns = F, cluster_row_slices = T, 
             column_labels = colnames(all)[1:14], show_row_names = FALSE, col = colmap, 
             row_split = all[,15], row_gap=unit(.05, "npc"), row_title_rot = 0)
draw(ht, padding = unit(c(2, 2, 25, 2), "mm"))

dev.off()


# Ion channel genes only

read.csv('/Users/admin/file.csv')