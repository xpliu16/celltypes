library(scrattch.mapping)
#install.packages('metap')
#library(metap)
library("readxl")

mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116"
#refFolder =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115/patchseq"
refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_116"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"
#patchseqFolder = r"(C:\Users\briank\Documents\Mapping\Macaque_MTGmapping\)"

#?????
#cell_list<-(read.csv2(r"(C:\Users\briank\Analysis\L5_cross_area\mapping to macaque MTG\cell_list_L5.csv)",header = TRUE,sep = ","))

#annoBG <- read_feather(file.path(refFolder, "anno.feather"))
#Expr.dat <- read_feather(file.path(refFolder, "data.feather"))
#AIT.anndata <- read_h5ad(file.path(refFolder,"NHP_BG_AIT115_taxonomy_entropy.h5ad"))
AIT.anndata <- loadTaxonomy(refFolder)
Expr.dat <- AIT.anndata$X
annoBG <- read_feather(file.path(refFolder, "anno.feather"))
#load(paste(refFolder,"/patchseq/QC_markers.RData",sep = ""))
#annoBG <- annoBG[match(Expr.dat$sample_id,annoBG$sample_id),]
annoBG <- annoBG[match(rownames(Expr.dat),annoBG$sample_id),]
#subclass<-c('L6 IT', 'L5 IT', 'L5 ET', 'L6 IT Car3', 'L6 CT', 'L5/6 NP', 'L6b',
#            'L4 IT')
#subclass<- unique(annoBG$Subclass_label[annoBG$level1.class_label!="NN"])      # BUT WE DON'T WANT NONNEURONAL CELL TYPES
# Striatal subclasses only (at least 5% of all cells are in dSTR or vSTR)
#subclass = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
#              "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
#              "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
#              "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
#              "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
#              "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
#              "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   
subclass = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome",
             "PVALB-COL19A1-ST18", "CCK-VIP-TAC3", "CCK-FBXL7",
             "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
             "D1-NUDAP", "TAC3-LHX8-PLPP4", "D1-ICj", "MEIS2", "SN_STH_GPe-MEIS2-OTX2", "SLC17A6")   

kpsubclass<-is.element(annoBG$Subclass_label,subclass)
annoBG<-annoBG[kpsubclass,]
Expr.dat<-Expr.dat[kpsubclass,]

# Color dictionary
colz <- annoBG$Subclass_color[match(subclass, annoBG$Subclass_label)]  
#col_table = data.frame(subclass, colz)
#write.csv(col_table, file.path(mappingFolder, 'colortable.csv'))

#patch_anno<-read_feather(paste(patchseqFolder, "anno2.feather",sep=""))
#load(file.path(mappingFolder,"NHP_BG_204_329_AIT115_ann_map_QC_sub.Rdata"))
#var<-load(file.path(mappingFolder,"NHP_BG_204_337_AIT115_ann_map_roi.Rdata"))
#load(file.path(mappingFolder,"NHP_BG_204_337_AIT115_ann_map_QC.Rdata"))
#annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_204_346_AIT115_ann_map_roi_QCpass.csv"))
#annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_AIT_116_RSC-204-370_sub_QC.csv"))
annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_AIT_116_RSC-204-378_roi_QC.csv"))

HANN_obj <- fromJSON(file.path(mappingFolder, "20240621_RSC-204-366_neurons_results.json"))
HANN_res <- do.call(cbind, HANN_obj$results)
orig_query <- read_h5ad(file.path(mappingFolder, '20240621_RSC-204-366_query.h5ad'))
HANN_cell_names <- orig_query$obs$exp_component_name
HANN_res <- cbind(HANN_cell_names, HANN_res)
rownames(HANN_res) = HANN_cell_names
HANN_res_sub <- HANN_res[annoNew$exp_component_name,]

patch_anno <- annoNew
#patch_anno <- merge(annoNew, HANN_res_sub, by.x='exp_component_name', by.y= 'row.names') 
#patch.dat<- read_feather(paste(patchseqFolder, "data.feather", sep =""))
load(paste0(data_dir, "/202411107_RSC-204-378_macaque_patchseq_star2.7_cpm.Rdata"))
patch.dat <- logCPM(cpmR) 

#patch_anno<- patch_anno[match(patch.dat$sample_id,patch_anno$sample_id),]
patch.dat <- patch.dat[,match(patch_anno$exp_component_name,colnames(patch.dat))] 
#kpSamp1 <- is.element(patch.dat$sample_id,cell_list[,1])
#patch_anno<-patch_anno[kpSamp1,]
#patch.dat<-patch.dat[kpSamp1,]
#patch_anno<-cbind(patch_anno,cell_list)
patch_anno$NMS_stringent<-patch_anno$marker_sum_norm_label > 0.5
patch.dat <- t(patch.dat)

#glial.genes <-as.vector(as.matrix(bind_cols(markers$Astro_on,markers$Endo_on,markers$Micro.PVM_on,markers$Oligo_on,markers$OPC_on,markers$VLMC_on,markers$OPC,
#                                            markers$VLMC)))
glial.genes <-as.vector(as.matrix(bind_cols(markers$Astrocytes_on,markers$Endothelial_on,markers$Microglia_on,markers$Oligos_on,markers$Oligos_Pre_on,
                                            markers$Ependyma_on)))
glial.exclude <- !(is.element(colnames(patch.dat),glial.genes))
patch.dat <- patch.dat[,glial.exclude]

sex.genes <- read.table(file.path(mappingFolder, "Jeremy_exclude_gene_lists/sex_genes.txt"),
                        sep = "")
sex.exclude <- !(is.element(colnames(patch.dat),sex.genes[[1]]))
patch.dat <- patch.dat[,sex.exclude]

genelist <- colnames(patch.dat)
genesSamp1 <- is.element(colnames(Expr.dat),genelist)
Expr.dat <- Expr.dat[,genesSamp1]

genelist2 <- colnames(Expr.dat)
genesSamp2 <- (is.element(colnames(patch.dat),genelist2))
patch.dat <- patch.dat[,genesSamp2]

# subset to patch dend genes
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
#dend <- readRDS(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
dend <- json_to_dend(fromJSON(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]]))
dendMarkers = unique(unlist(get_dend_markers(dend)))
genesSamp1 <- is.element(colnames(Expr.dat),dendMarkers)
Expr.dat <- Expr.dat[,genesSamp1]
genesSamp2 <- (is.element(colnames(patch.dat),dendMarkers))
patch.dat <- patch.dat[,genesSamp2]

#Expr.dat<-as.matrix(Expr.dat[,names(Expr.dat)!="sample_id"])
#Expr.dat<-as.matrix(Expr.dat[,names(Expr.dat)!="cell_id"])   # These are non unique?
#rownames(Expr.dat)<-annoBG$sample_id
Expr.dat<-t(Expr.dat)
#patch.dat<-as.matrix(patch.dat[,names(patch.dat)!="cell_id"])

rownames(patch.dat)<-patch_anno$exp_component_name
patch.dat<-t(patch.dat)    # BOTH INTO GENE X CELL
#patch.dat<-log2(patch.dat+1)    # CHECK MAX(PATCH.DAT) TO SEE IF YOU NEED TO LOG TRANSFORM

Expr.dat<-Expr.dat[rownames(patch.dat),colnames(Expr.dat)]  # Indexing by names - why not just ,] ?
platformGn = abs(rowMeans(Expr.dat)-rowMeans(patch.dat))>=3
expressedGn = (rowSums(Expr.dat>=1)>(0.01*dim(Expr.dat)[2]))&(rowSums(patch.dat>=1)>(0.01*dim(patch.dat)[2]))
keepGenes <- (!(platformGn))&(expressedGn)
mean(keepGenes)   # Percentage we're keeping

dataBG_all<-Expr.dat
annoBG_all<-annoBG
dataBG_all_PS<-patch.dat
annoBG_all_PS<-patch_anno

annoBG_all_PS$Prep = NA
annoBG_all_PS$Prep[annoBG_all_PS$cell_specimen_project == "qIVSCC-METa"] = "acute"
annoBG_all_PS$Prep[annoBG_all_PS$cell_specimen_project == "qIVSCC-METc"] = "cultured"

link <- read_excel(file.path(mappingFolder, 'HMBA_BG_consensus_annotation.xlsx'), sheet = 'linking_table')
link_vals <- link$level4_group
keys <- link$previous_nomenclature
names(link_vals) <- keys

annoBG_all$Subclass_label = link_vals[annoBG_all$Subclass_label]
annoBG_all_PS$Subclass_Corr_former = annoBG_all_PS$Subclass_Corr
annoBG_all_PS$Subclass_Corr = link_vals[annoBG_all_PS$Subclass_Corr]
annoBG_all_PS$Subclass_Tree_former = annoBG_all_PS$Subclass_Tree
annoBG_all_PS$Subclass_Tree = link_vals[annoBG_all_PS$Subclass_Tree]

## Basic data and metadata set-up

brain.data     <- cbind(dataBG_all[keepGenes,],dataBG_all_PS[keepGenes,])  # Include only genes subsetted above
brain.metadata <- data.frame(set=c(rep("Taxonomy",dim(dataBG_all)[2]),rep("Patch-seq",dim(dataBG_all_PS)[2])),
#                             subclass = c(annoBG_all$Subclass_label,annoBG_all_PS$Subclass.assignment),   # Temporarily changed to Corr
                             subclass = c(annoBG_all$Subclass_label,annoBG_all_PS$Subclass_Corr), 
                             QC = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$Norm_Marker_Sum.0.4_label), 
                             QC2 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$NMS_stringent),
                             QC3 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$library_prep_pass_fail),
                             QC4 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$percent_reads_aligned_to_introns>25),
                             QC5 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$score.Corr>0.55),
                             QC6 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$Genes.Detected<8000),
                             QC6 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$amplified_quantity_ng),
#                             subclass_color = c(annoBG_all$Subclass_color, annoBG_all$Subclass_color[match(annoBG_all_PS$Subclass.assignment, annoBG_all$Subclass_label)]),   # TEMPORARILY CHANGED TO CORR
                             subclass_color = c(annoBG_all$Subclass_color, annoBG_all$Subclass_color[match(annoBG_all_PS$Subclass_Corr, annoBG_all$Subclass_label)]), 
                             area =c(rep("BG",dim(dataBG_all)[2]),annoBG_all_PS$Region),
                             subclass_tree =c(annoBG_all$Subclass_label,annoBG_all_PS$Subclass_Tree),
                             subclass_color_tree = c(annoBG_all$Subclass_color, annoBG_all$Subclass_color[match(annoBG_all_PS$Subclass_Tree, annoBG_all$Subclass_label)]),
                             prep = c(rep("acute",dim(dataBG_all)[2]),annoBG_all_PS$Prep),
                             virus = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$Virus),
                             virus4609 = c(rep("none",dim(dataBG_all)[2]), (annoBG_all_PS$Virus=="CN4609") & (annoBG_all_PS$creCell=="Positive")),
                             virus5050 = c(rep("none",dim(dataBG_all)[2]), (annoBG_all_PS$Virus=="CN5050") & (annoBG_all_PS$creCell=="Positive")),                             
                             creCell = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$creCell))
rownames(brain.metadata) <- colnames(brain.data)

## Construct data set lists
brain      <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)
brain.list <- SplitObject(object = brain, split.by = "set")

# beta score which indicates the binaryness of a gene across clusters. 
# High scores (near 1) indicate that a gene is either on or off in nearly all cells of every cluster. 
cl<- setNames(annoBG_all$cluster_label,colnames(dataBG_all))
propExpr  <- get_cl_prop(dataBG_all[keepGenes,],cl)
betaScore <- getBetaScore(propExpr)
betaOut   <- data.frame(Gene=rownames(dataBG_all)[keepGenes],BetaScore=betaScore)
betaOut   <- betaOut[order(-betaScore),]

nGenes      <- min(2000, length(betaOut$Gene))
varFeatures <- betaOut$Gene[1:nGenes]

for (i in 1:length(x = brain.list)) {
  VariableFeatures(brain.list[[i]]) <- varFeatures
}

dims = 30
#brain.list <- lapply(X = brain.list, FUN = function(x) {
# x <- NormalizeData(x)
#x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})

# select features that are repeatedly variable across datasets for integration
#features <- SelectIntegrationFeatures(object.list = brain.list)
brain.anchors <- FindIntegrationAnchors(object.list = brain.list, dims =1:dims,
                                        anchor.features = varFeatures,verbose =FALSE,k.anchor = 5, k.filter = 100)
# this command creates an 'integrated' data assay
brain.combined <- IntegrateData(anchorset = brain.anchors,dims=1:dims,verbose = FALSE)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(brain.combined) <- "integrated"

# Run the standard workflow for visualization and clustering

brain.combined <- ScaleData(brain.combined, verbose = FALSE)
brain.combined <- RunPCA(brain.combined, npcs = dims, verbose = FALSE)

#brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:dims)
brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:dims, 
                          n.neighbors = 180, n.components = 2, min.dist = 0.3)
brain.combined <- FindNeighbors(brain.combined, reduction = "pca", dims = 1:dims)
brain.combined <- FindClusters(brain.combined, resolution = 1)

# Visualization

xl <- range(brain.combined@reductions$umap@cell.embeddings[,1])
yl <- range(brain.combined@reductions$umap@cell.embeddings[,2])

mappingFolder<- file.path(mappingFolder,"Int_UMAP_378_STR_roi")
dir.create(mappingFolder, showWarnings = FALSE)
#dir.create(mappingFolder, showWarnings=FALSE)

jpeg(file.path(mappingFolder,'204_378_integrated_umap_dendGn_sub_180neigh.jpg'), width = 2000, height = 1200, quality = 100)

DimPlot(brain.combined, reduction = "umap", split.by = "set",pt.size = 1)
umap = as.matrix(brain.combined@reductions[["umap"]]@cell.embeddings)
umap1 = (umap[,1] > 0)&(umap[,1]<3)
umap = umap[umap1,]
umap2 = (umap[,2] > -6)&(umap[,2]<-4.5)
umap = umap[umap2,]

row.names(umap)

dev.off()

# For all subclasses
#jpeg(file.path(mappingFolder,'204_337_integrated_umap2_dendGn.jpg'), width = 2500, height = 1600, quality = 100)
# For STR subclasses only
jpeg(file.path(mappingFolder,'204_378_integrated_umap2_dendGn_sub.jpg'), width = 2000, height = 1600, quality = 100)

p1 <- DimPlot(brain.combined, reduction = "umap", group.by = "set", label.size = 1, pt.size = 1, cols = c("red", "grey"), raster = F)+xlim(xl) + ylim(yl)
p2 <- DimPlot(brain.combined, reduction = "umap", group.by = "subclass", cells=colnames(dataBG_all),pt.size = .5, raster=F)+ggtitle("FACS")+xlim(xl) + ylim(yl) +guides(color = guide_legend(override.aes = list(size=4), ncol=1))
p3 <- DimPlot(brain.combined, group.by = "subclass",  reduction = "umap",
              pt.size = 1.5, label=FALSE, label.size = 2,cells=colnames(dataBG_all_PS), raster=F)+ ggtitle("Patch-seq (Corr)")+xlim(xl) + ylim(yl)
p4 <- DimPlot(brain.combined, group.by = "area",  reduction = "umap",
              pt.size = 1.5, label=FALSE, label.size = 2,cells=colnames(dataBG_all_PS), raster=F)+ ggtitle("ROI Area")+xlim(xl) + ylim(yl)
p5 <- DimPlot(brain.combined, group.by = "subclass_tree",  reduction = "umap",
              pt.size = 1.5, label=FALSE, label.size = 2,cells=colnames(dataBG_all_PS), raster=F)+ ggtitle("Patch-seq (Tree)")+xlim(xl) + ylim(yl)
p6 <- DimPlot(brain.combined, group.by = "prep",  reduction = "umap",
              pt.size = 1.5, label=FALSE, label.size = 2,cells=colnames(dataBG_all_PS), raster=F)+ ggtitle("Prep")+xlim(xl) + ylim(yl)
plot_grid(p1, p2, p3,p4,p5,p6, ncol=2)

dev.off()


#unique(brain.combined$subclass_corr[brain.combined$set=="PatchSeq"])
#unique(brain.combined$subclass[brain.combined$set=="PatchSeq"])
#cols = c('V0' = 'red', 'V6' = 'grey', 'V8' = 'grey')

jpeg(file.path(mappingFolder,'204_378_integrated_umap3_dendGn_sub.jpg'), width = 2400, height = 1200,
     pointsize = 12, quality = 100)
plot_grid(p1,p2)
p1+p2    # WHAT DOES THIS DO?
dev.off()

#install.packages("pals")
#library(pals)
colz <- DiscretePalette(n = length(unique(annoBG_all$Subclass_label)), palette = "polychrome")
jpeg(file.path(mappingFolder,'204_378_integrated_umap4_dendGn_lab_repel_sub.jpg'), width = 2600, height = 1200,
     pointsize = 12, quality = 100)
p2 <- DimPlot(brain.combined, reduction = "umap", group.by = "subclass", cells=colnames(dataBG_all),pt.size = .5, cols=colz, label=T, repel=T, raster=F)+ggtitle("FACS")+xlim(xl) + ylim(yl)
#LabelClusters(plot = p2, id = 'subclass')
plot_grid(p1,p2)
p1+p2    # WHAT DOES THIS DO?
dev.off()

# Used for figures

sc_names = sort(unique(brain.combined@meta.data$subclass))
#colz <- brain.combined@meta.data$subclass_color_tree[match(sc_names, brain.combined@meta.data$subclass)]  
colz <- DiscretePalette(n = length(unique(annoBG_all$Subclass_label)), palette = "polychrome", shuffle=TRUE)
#jpeg(file.path(mappingFolder,'204_359_integrated_umap4_dendGn_nolab2_sub.jpg'), width = 3500, height = 1500,
#     pointsize = 12, quality = 100)
png(file.path(mappingFolder,'204_378_integrated_umap4_dendGn_nolab2_sub.png'), width = 3100, height = 1500, pointsize = 12)
p2 <- DimPlot(brain.combined, reduction = "umap", group.by = "subclass", cells=colnames(dataBG_all),pt.size = 1, cols=colz, raster=F)+
  xlim(xl) + ylim(yl) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1) ) 
#LabelClusters(plot = p2, id = 'subclass')
p1<- p1 + theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45), axis.title=element_text(size=45, face="bold", margin=0.1), legend.position="none")
p2<- p2 + theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45), axis.title=element_text(size=45, face="bold", margin=0.1), legend.text = element_text(size=25))
plot_grid(p1,p2)
p1+p2    # WHAT DOES THIS DO?
dev.off()

#umap_tx = data.frame(umap1 = umap1, umap2 = umap2, colz = colz)
#jpeg(file.path(mappingFolder,'204_337_integrated_umap4_dendGn_nolab_sub_ggplot.jpg'), width = 2400, height = 1200,
#     pointsize = 12, quality = 100)
#ggplot(umap_tx, aes(x=umap1, y=umap2, color=colz)) + geom_point() 
#dev.off()

umap1 <- brain.combined@reductions$umap@cell.embeddings[,1]
umap2 <- brain.combined@reductions$umap@cell.embeddings[,2]
brain.metadata <- cbind(brain.metadata, umap1, umap2)

main="UMAP Visualization of Corr mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.6  # 0.4
pch=19
cex.main=1
cex.legend=1

inds_ps = brain.metadata$set == "Patch-seq"
inds_tx = brain.metadata$set == "Taxonomy"

jpeg(file.path(mappingFolder,'204_378_integrated_umap_Corr_Tree_dendGn_sub.jpg'), width = 1500, height = 600,
     pointsize = 12, quality = 100)
par(mfrow=c(1,3))
plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col=brain.metadata$subclass_color[inds_tx], cex=cex, pch=pch)
mtext(side=3, "BG Taxonomy", cex=cex.main)
labels.u <- sort(unique(brain.metadata$subclass))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- brain.metadata$subclass_color[match(labels.u,brain.metadata$subclass)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col="grey84", cex=cex, pch=pch)
points(umap1[inds_ps], umap2[inds_ps], col=brain.metadata$subclass_color[inds_ps], cex=cex, pch=pch)
mtext(side=3, "Patchseq Cells (Corr)", cex=cex.main)

plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col="grey84", cex=cex, pch=pch)
points(umap1[inds_ps], umap2[inds_ps], col=brain.metadata$subclass_color_tree[inds_ps], cex=cex, pch=pch)
mtext(side=3, "Patchseq Cells (Tree)", cex=cex.main)

dev.off()

jpeg(file.path(mappingFolder,'204_378_integrated_umap_QC_dendGn_sub.jpg'), width = 2400, height = 1600, quality = 100)
p7 <- DimPlot(brain.combined, reduction = "umap", pt.size = 1, cells=colnames(dataBG_all_PS), group.by = "QC") + ggtitle("NMS > 0.4")
#p8 <- DimPlot(brain.combined, reduction = "umap", cells=colnames(dataBG_all_PS), group.by = "QC2") + ggtitle("NMS stringent")
p8 <- DimPlot(brain.combined, reduction = "umap",pt.size = 1, cells=colnames(dataBG_all_PS), group.by = "QC5") + ggtitle("score.Corr > 0.55")
p9 <- DimPlot(brain.combined, reduction = "umap", pt.size = 1, cells=colnames(dataBG_all_PS), group.by = "QC4") + ggtitle("Percent Reads Intronic > 25%")

plot_grid(p1, p2, p6, p7, p8, p9, nrow=2)
dev.off()

jpeg(file.path(mappingFolder,'204_378_integrated_umap_QC2_dendGn_sub.jpg'), width = 2400, height = 1600, quality = 100)
p10 <- DimPlot(brain.combined, reduction = "umap", pt.size = 1, cells=colnames(dataBG_all_PS), group.by = "QC3") + ggtitle("Library_prep_pass_fail")
p11 <- DimPlot(brain.combined, reduction = "umap", pt.size = 1, cells=colnames(dataBG_all_PS), group.by = "QC6") + ggtitle("Genes detected < 8000")
p12 <- DimPlot(brain.combined, reduction = "umap", pt.size = 1, cells=colnames(dataBG_all_PS), group.by = "QC2") + ggtitle("NMS_stringent")
plot_grid(p1, p2, p7, p10, p11, p12, nrow=2)
dev.off()

# By virus
jpeg(file.path(mappingFolder,'204_378_integrated_umap_virus_roi.jpg'), width = 2400, height = 1000, quality = 100)
p7 <- DimPlot(brain.combined, group.by = "virus",  reduction = "umap",
              pt.size = 1.5, label=FALSE, label.size = 2,cells=colnames(dataBG_all_PS), raster=F)+ ggtitle("Prep")+xlim(xl) + ylim(yl)
dev.off()

jpeg(file.path(mappingFolder,'204_378_integrated_umap_virus4609_roi_b.jpg'), width = 2000, height = 1000, quality = 100)

par(mfrow=c(1,2))
plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col=brain.metadata$subclass_color[inds_tx], cex=cex, pch=pch)
mtext(side=3, "BG Taxonomy", cex=cex.main)
labels.u <- sort(unique(brain.metadata$subclass))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- brain.metadata$subclass_color[match(labels.u,brain.metadata$subclass)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col="grey84", cex=cex, pch=pch)
points(umap1[inds_ps], umap2[inds_ps], col="cyan", cex=cex, pch=pch)
inds_ps_virus = brain.metadata$virus4609 == TRUE
points(umap1[inds_ps_virus], umap2[inds_ps_virus], col= "red", cex=cex, pch=pch)
mtext(side=3, "Patchseq Cells (CN4609)", cex=cex.main)

dev.off()

jpeg(file.path(mappingFolder,'204_378_integrated_umap_virus5050_roi_b.jpg'), width = 2400, height = 1000, quality = 100)
par(mfrow=c(1,2))
plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col=brain.metadata$subclass_color[inds_tx], cex=cex, pch=pch)
mtext(side=3, "BG Taxonomy", cex=cex.main)
labels.u <- sort(unique(brain.metadata$subclass))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- brain.metadata$subclass_color[match(labels.u,brain.metadata$subclass)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

plot(xl, yl, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xl[1], yl[1], xl[2], yl[2], border="#aaaaaa", lwd=0.25) 
points(umap1[inds_tx], umap2[inds_tx], col="grey84", cex=cex, pch=pch)
points(umap1[inds_ps], umap2[inds_ps], col="cyan", cex=cex, pch=pch)
inds_ps_virus = brain.metadata$virus5050 == TRUE
points(umap1[inds_ps_virus], umap2[inds_ps_virus], col= "red", cex=cex, pch=pch)
mtext(side=3, "Patchseq Cells (CN5050)", cex=cex.main)

dev.off()

# High genes detected
brain.metadata_sub <- brain.metadata[inds_ps,]
unique(brain.metadata_sub$subclass[!(as.numeric(brain.metadata_sub$QC6=="TRUE"))])


DimPlot(brain.combined, reduction = "umap", split.by = "set",pt.size = 1, group.by = "subclass")

DimPlot(brain.combined, reduction = "umap", split.by = "set",pt.size = 1, group.by = "QC")
DimPlot(brain.combined, reduction = "umap", split.by = "set",pt.size = 1, group.by = "QC2")

DefaultAssay(brain.combined) <- "RNA"
nk.markers <- FindConservedMarkers(brain.combined, ident.1 = "12", grouping.var = "set", verbose = FALSE)
head(nk.markers)

ETmarker<-as.vector(markers$L5.ET_on)
FeaturePlot(brain.combined, features = c("GRM7","TOX2","STXBP6","TRPC4","ERG","SNCG"), min.cutoff = "q25")


###
#Grab all PVALBs and analyze split

PV = brain.metadata[brain.metadata$subclass_corr == "PVALB-COL19A1-ST18",]
xl_PV = umap1[brain.metadata$subclass_corr == "PVALB-COL19A1-ST18"]
yl_PV = umap2[brain.metadata$subclass_corr == "PVALB-COL19A1-ST18"]
inds = PV$set == "Patch-seq"
PV = PV[inds,]
xl_PV = xl_PV[inds]
yl_PV = yl_PV[inds]

jpeg(file.path(mappingFolder,'PVALB_xl_hist.jpg'), width = 500, height = 350, quality = 100)
hist(xl_PV, 20)
dev.off()

# Cut at -1.5

PV_l = PV[xl_PV<(-1.5),]
PV_r = PV[xl_PV>=(-1.5),]
PV_l_names <- rownames(PV_l)
PV_r_names <- rownames(PV_r)

table(patch_anno[patch_anno$exp_component_name %in% PV_l_names,'Cluster_Corr'])
#21_IN 22_IN 24_IN 25_IN 27_IN 28_IN 
#10    14     2     1     3     3

table(patch_anno[patch_anno$exp_component_name %in% PV_r_names,'Cluster_Corr'])
#21_IN 22_IN 23_IN 24_IN 25_IN 27_IN 28_IN 
#25     1     2     3     1    18     2 

table(patch_anno[patch_anno$exp_component_name %in% PV_l_names,'Region'])
#STRd_Ca  STRd_Pu   STRdCa   STRdPu STRv_ACB  STRvACB 
#15        7        4        3        3        1 

table(patch_anno[patch_anno$exp_component_name %in% PV_r_names,'Region'])
#STRd_Ca  STRd_Pu   STRdCa   STRdPu STRv_ACB  STRvACB 
#20       13        2        5        6        6 

table(patch_anno[patch_anno$exp_component_name %in% PV_l_names,'Region'])
table(patch_anno[patch_anno$exp_component_name %in% PV_l_names,'Region'])


# Taxonomy cells
PV = brain.metadata[brain.metadata$subclass_corr == "PVALB-COL19A1-ST18",]
xl_PV = umap1[brain.metadata$subclass_corr == "PVALB-COL19A1-ST18"]
yl_PV = umap2[brain.metadata$subclass_corr == "PVALB-COL19A1-ST18"]
inds = PV$set == "Taxonomy"
PV = PV[inds,]
xl_PV = xl_PV[inds]
yl_PV = yl_PV[inds]

PV_l_tax = PV[xl_PV<(-1.0),]
PV_r_tax = PV[xl_PV>=(-1.0),]
PV_l_tax_names <- rownames(PV_l_tax)
PV_r_tax_names <- rownames(PV_r_tax)

jpeg(file.path(mappingFolder,'PVALB_xl_hist_tax.jpg'), width = 500, height = 350, quality = 100)
hist(xl_PV, 20)
dev.off()

table(annoBG[annoBG$sample_id %in% PV_l_tax_names,'Cluster_label'])
# 21_IN 22_IN 23_IN 24_IN 25_IN 26_IN 27_IN 28_IN 
#  505   603   181   302   115   599   481   403 
table(annoBG[annoBG$sample_id %in% PV_r_tax_names,'Cluster_label'])
#21_IN 22_IN 23_IN 24_IN 25_IN 26_IN 27_IN 28_IN 
# 495   397   161   226   103   401   519   161 

table(annoBG[annoBG$sample_id %in% PV_l_tax_names,'Region_label'])
#Macaque CaB  Macaque CaH  Macaque CaT  Macaque GPe  Macaque GPi  Macaque NAC 
#402         1046          423           26           20          451 
#Macaque PuC Macaque PuPV  Macaque PuR 
#304          272          245 

table(annoBG[annoBG$sample_id %in% PV_r_tax_names,'Region_label'])
Region_label
#Macaque CaB  Macaque CaH  Macaque CaT  Macaque GPe  Macaque NAC  Macaque PuC 
#500          534          239           65          459          195 
#Macaque PuPV  Macaque PuR  Macaque STH 
#156          314            1 

table(annoBG[annoBG$sample_id %in% PV_l_tax_names,'Species_label'])
# Macaca mulatta Macaca nemestrina 
#     1352              1837 


save(PV_l_names, PV_r_names, PV_l_tax_names, PV_r_tax_names, file=file.path(mappingFolder, 'PVALB_left_right.Rdata'))

load(file=file.path(mappingFolder, 'Int_UMAP/PVALB_left_right.Rdata'))
#data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"
#load(paste0(data_dir, "/20240520_RSC-204-370_macaque_patchseq_star2.7_cpm.Rdata"))
#patch.dat <- logCPM(cpmR) 
annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_AIT_116_RSC-204-370_sub_QC.csv"))
patch_anno <- annoNew

table(patch_anno[patch_anno$exp_component_name %in% PV_l_names,'gender'])
table(patch_anno[patch_anno$exp_component_name %in% PV_r_names,'gender'])
