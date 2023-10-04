#refFolder <- "/home/xiaoping.liu/scrattch/reference/NHP_BG_AIT_114"
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/"
#dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

library(umap)
library(scrattch.mapping)

#load(file.path(mappingFolder,"NHP_BG_RSC_204_324_mapping.Rdata"))
#load(file.path(mappingFolder,"NHP_BG_RSC_204_324_ann_map.Rdata"))
load(file.path(mappingFolder,"NHP_BG_204_337_AIT115_mapping2.Rdata"))
#load(file.path(mappingFolder,"NHP_BG_204_337_AIT115ann_map_full.Rdata"))    # annotations_mapped with all samples

load(file.path(mappingFolder,"NHP_BG_204_337_AIT115_ann_map_QC_full.Rdata"))  # anno_new with basic QC'ed samples

load(paste0(data_dir, "/20230807_RSC-204-337_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20230807_RSC-204-337_macaque_patchseq_star2.7_samp.dat.Rdata"))

query.metadata <- samp.dat
counts      <- cpmR   # Genes are rows, samples are columns

query.counts   <- counts
query.data   <- logCPM(query.counts)

# Put annotations and counts in the same order
query.metadata <- query.metadata[match(colnames(query.data),query.metadata$exp_component_name),] 
rownames(query.metadata) <- query.metadata$exp_component_name  

#Hack
rownames(query.mapping) = rownames(query.metadata)

#AIT.anndata <- read_h5ad(file.path(refFolder,"AIT_114_taxonomy.h5ad"))
AIT.anndata <- loadTaxonomy(refFolder)
# Assign levels
clusters <- unique(AIT.anndata$uns$clusterInfo$cluster_label)
subclass_lavels <- unique(AIT.anndata$uns$clusterInfo$level3.subclass_label)

annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = 0, all=TRUE) 
annotations_mapped <- annotations_mapped[match(rownames(query.metadata),annotations_mapped$exp_component_name),] 
rownames(annotations_mapped) <- annotations_mapped$exp_component_name 

#annotations_mapped$cluster <- factor(query.mapping$cluster, levels=clusters)  # Make into discrete levels
#inds1 = grepl("STR",annotations_mapped$roi)
inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annotations_mapped$roi), TRUE,FALSE)
inds2 = annotations_mapped$library_prep_pass_fail == "Pass"
inds3 = annotations_mapped$Genes.Detected >= 1000
inds4 = annotations_mapped$percent_reads_aligned_total >= 50
anno_mapped_sub = annotations_mapped[inds1&inds2&inds3&inds4,]  # JUST USE INDS1 HERE TO MATCH UMAP OF 946 samples, then apply basic QC below
query.data_sub = query.data[,inds1&inds2&inds3&inds4]
# check n cells

#query.counts_t <- t(query.counts)
#mapping.data <- query.counts_t
query.data_t <- t(query.data_sub)
mapping.data <- query.data_t
mapping.labels_Corr <- anno_mapped_sub$level3.subclass_Corr
mapping.labels_Tree <- anno_mapped_sub$level3.subclass_Tree
# check rows are aligned - yes

#mapping.umap <- umap(mapping.data)
#save(mapping.umap, file=file.path(mappingFolder,"NHP_BG_204_337_AIT115_umap_roi.Rdata"))
load(file=file.path(mappingFolder,"NHP_BG_204_337_AIT115_umap_roi.Rdata"))

# Optional: Further subset by score.Corr
inds5 <- anno_mapped_sub$score.Corr >= 0.6
mapping.umap$layout <- mapping.umap$layout[inds5,]
anno_mapped_sub <- anno_mapped_sub[inds5,]
mapping.labels_Corr <- mapping.labels_Corr[inds5]
mapping.labels_Tree <- mapping.labels_Tree[inds5]

# Optional: apply basic QC
inds3 = anno_mapped_sub$Genes.Detected >= 1000
inds4 = anno_mapped_sub$percent_reads_aligned_total >= 50
mapping.umap$layout <- mapping.umap$layout[inds3&inds4,]
anno_mapped_sub <- anno_mapped_sub[inds3&inds4,]
mapping.labels_Corr <- mapping.labels_Corr[inds3&inds4]
mapping.labels_Tree <- mapping.labels_Tree[inds3&inds4]

mapping.colors_Tree <- 
  AIT.anndata$uns$clusterInfo$level3.subclass_color[match(mapping.labels_Tree, AIT.anndata$uns$clusterInfo$level3.subclass_label)]

mapping.colors_Corr <- 
  AIT.anndata$uns$clusterInfo$level3.subclass_color[match(mapping.labels_Corr, AIT.anndata$uns$clusterInfo$level3.subclass_label)]

layout <- mapping.umap$layout

#jpeg(file.path(mappingFolder,'204_337_umap_Tree_roi_nGenes_percAligned.jpg'), quality = 100)
jpeg(file.path(mappingFolder,'204_337_umap_Tree_roi.jpg'), quality = 100)
main="UMAP Visualization of Tree mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

points(layout[,1], layout[,2], col=mapping.colors_Tree, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_Tree))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_Tree[match(labels.u,mapping.labels_Tree)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()
# save figure


jpeg(file.path(mappingFolder,'204_324_umap_Corr_scoreQC.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_Corr))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_Corr[match(labels.u,mapping.labels_Corr)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_324_umap_Corr_MSN.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr_MSN <- mapping.colors_Corr
mapping.colors_Corr_MSN [anno_mapped_sub$class_Corr != "MSN"] = "#808080"

points(layout[,1], layout[,2], col=mapping.colors_Corr_MSN, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_Corr[anno_mapped_sub$class_Corr == "MSN"]))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_Corr_MSN[match(labels.u,mapping.labels_Corr)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_324_umap_Tree_MSN.jpg'), quality = 100)
main="UMAP Visualization of Tree mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Tree_MSN <- mapping.colors_Tree
mapping.colors_Tree_MSN [anno_mapped_sub$class_Tree != "MSN"] = "#808080"

points(layout[,1], layout[,2], col=mapping.colors_Tree_MSN, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_Tree[anno_mapped_sub$class_Tree == "MSN"]))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_Tree_MSN[match(labels.u,mapping.labels_Tree)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

mapping.labels_class_Tree <- anno_mapped_sub$class_Tree
mapping.colors_class_Tree <- 
  AIT.anndata$uns$clusterInfo$class_color[match(mapping.labels_class_Tree, AIT.anndata$uns$clusterInfo$class_label)]

jpeg(file.path(mappingFolder,'204_324_umap_Tree_class.jpg'), quality = 100)
main="UMAP Visualization of Tree mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

points(layout[,1], layout[,2], col=mapping.colors_class_Tree, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_class_Tree))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_class_Tree[match(labels.u,mapping.labels_class_Tree)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

mapping.labels_class_Corr <- anno_mapped_sub$class_Corr
mapping.colors_class_Corr <- 
  AIT.anndata$uns$clusterInfo$class_color[match(mapping.labels_class_Corr, AIT.anndata$uns$clusterInfo$class_label)]

jpeg(file.path(mappingFolder,'204_324_umap_Corr_class.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

points(layout[,1], layout[,2], col=mapping.colors_class_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_class_Corr))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_class_Corr[match(labels.u,mapping.labels_class_Corr)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


jpeg(file.path(mappingFolder,'204_324_umap_Corr_MSN_hybrid.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr_MSN <- mapping.colors_Corr
mapping.colors_Corr_MSN [anno_mapped_sub$class_Corr != "MSN"] = "#808080"
mapping.colors_Corr_MSN [anno_mapped_sub$subclass_Corr == "D1D2 Hybrid"] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr_MSN, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_Corr[anno_mapped_sub$class_Corr == "MSN"]))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_Corr_MSN[match(labels.u,mapping.labels_Corr)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_324_umap_Tree_MSN_hybrid.jpg'), quality = 100)
main="UMAP Visualization of Tree mapped NHP_BG Patch-seq cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Tree_MSN <- mapping.colors_Tree
mapping.colors_Tree_MSN [anno_mapped_sub$class_Tree != "MSN"] = "#808080"
mapping.colors_Tree_MSN [anno_mapped_sub$subclass_Tree == "D1D2 Hybrid"] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Tree_MSN, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- sort(unique(mapping.labels_Tree[anno_mapped_sub$class_Tree == "MSN"]))
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- mapping.colors_Tree_MSN[match(labels.u,mapping.labels_Tree)]
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()



jpeg(file.path(mappingFolder,'204_337_AIT115_umap_library_prep.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: library_prep_pass_fail"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$library_prep_pass_fail == "Pass"] = "#808080"
mapping.colors_Corr [anno_mapped_sub$library_prep_pass_fail == "Fail"] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("Pass", "Fail")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_percent_intronic_15.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: percent_reads_aligned_to_introns"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$percent_reads_aligned_to_introns > 15] = "#808080"
mapping.colors_Corr [anno_mapped_sub$percent_reads_aligned_to_introns <= 15] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("> 15%", "<= 15%")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


jpeg(file.path(mappingFolder,'204_337_AIT115_umap_scoreCorr2.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: score.Corr (mapping confidence)"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$score.Corr > 0.4] = "#808080"
mapping.colors_Corr [anno_mapped_sub$score.Corr <= 0.4] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("> 0.4", "<= 0.4")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'scoreCorr_dist.jpg'))
hist_info <- hist(anno_mapped_sub$score.Corr, 
                  main = "score.Corr distribution",
                  freq = TRUE, plot = TRUE)
dev.off()


jpeg(file.path(mappingFolder,'204_337_AIT115_umap_nGenes.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Number of genes detected"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$Genes.Detected < 8000] = "#808080"
mapping.colors_Corr [anno_mapped_sub$Genes.Detected >= 8000] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("< 8000", ">= 8000")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_nGenes_1000.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Number of genes detected"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$Genes.Detected > 1000] = "#808080"
mapping.colors_Corr [anno_mapped_sub$Genes.Detected <= 1000] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("> 1000", "<= 1000")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


jpeg(file.path(mappingFolder,'204_337_AIT115_umap_400bp.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Percent > 400bp"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$percent_cdna_longer_than_400bp > 0.3] = "#808080"
mapping.colors_Corr [anno_mapped_sub$percent_cdna_longer_than_400bp <= 0.3] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("> 0.3", "<= 0.3")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


jpeg(file.path(mappingFolder,'204_337_AIT115_umap_amplified_quantity.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Amplified quantity (ng)"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
inds = (anno_mapped_sub$amplified_quantity_ng > 4) & (anno_mapped_sub$amplified_quantity_ng < 90)
mapping.colors_Corr [inds] = "#808080"
mapping.colors_Corr [!inds] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("4 < x < 90", "others")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_amplified_quantity_4.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Amplified quantity (ng)"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
inds = (anno_mapped_sub$amplified_quantity_ng > 4)
mapping.colors_Corr [inds] = "#808080"
mapping.colors_Corr [!inds] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("Amplified quantity > 4 ng", "Amplified quantity <= 4 ng")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_amplified_quantity_90.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Amplified quantity (ng)"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
inds = (anno_mapped_sub$amplified_quantity_ng < 90)
mapping.colors_Corr [inds] = "#808080"
mapping.colors_Corr [!inds] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("Amplified quantity < 90 ng", "Amplified quantity >= 90 ng")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

quality_score <- annoNew[match(anno_mapped_sub$exp_component_name, annoNew$exp_component_name), 
               "quality_score_label"] 
jpeg(file.path(mappingFolder,'204_337_AIT115_umap_quality_score.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Quality Score (Pavlidis lab)"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
inds = (quality_score > 0.2)
mapping.colors_Corr [inds] = "#808080"
mapping.colors_Corr [!inds] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("Quality Score > 0.2", "Quality Score <= 0.2")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


jpeg(file.path(mappingFolder,'204_337_AIT115_umap_culture_acute.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Acute vs. Cultured"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$cell_specimen_project == 'qIVSCC-METa'] = "#808080"
mapping.colors_Corr [anno_mapped_sub$cell_specimen_project == 'qIVSCC-METc'] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("Acute", "Cultured")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

NMS <- annoNew[match(anno_mapped_sub$exp_component_name, annoNew$exp_component_name), 
               "Norm_Marker_Sum.0.4_label"] 
jpeg(file.path(mappingFolder,'204_337_AIT115_umap_NMS.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Normalized Marker Sum"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [NMS == 'TRUE'] = "#808080"
mapping.colors_Corr [NMS == 'FALSE'] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("NMS > 0.4", "NMS <= 0.4")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

NMS <- annoNew[match(anno_mapped_sub$exp_component_name, annoNew$exp_component_name), 
               "marker_sum_norm_label"] 
NMS <- NMS > 0.6

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_NMS06.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Normalized Marker Sum"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [NMS == 'TRUE'] = "#808080"
mapping.colors_Corr [NMS == 'FALSE'] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("NMS > 0.6", "NMS <= 0.6")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


contamsum <- annoNew[match(anno_mapped_sub$exp_component_name, annoNew$exp_component_name), 
                     "contam_sum_label"] 
jpeg(file.path(mappingFolder,'contam_sum_dist.jpg'))
hist_info <- hist(contamsum, 
                  main = "Contam_sum",
                  freq = TRUE, plot = TRUE)
dev.off()

contamsum <- contamsum < 4.0

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_contamsum.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: Contam_sum"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [contamsum == 'TRUE'] = "#808080"
mapping.colors_Corr [contamsum == 'FALSE'] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("contam_sum < 4.0", "contam_sum >= 4.0")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


jpeg(file.path(mappingFolder,'204_337_AIT115_umap_percaligned25.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: percent_reads_aligned"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$percent_reads_aligned_total > 25] = "#808080"
mapping.colors_Corr [anno_mapped_sub$percent_reads_aligned_total <= 25] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("percent_reads_aligned_total > 25", "percent_reads_aligned_total <= 25")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_percaligned50.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells: percent_reads_aligned"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [anno_mapped_sub$percent_reads_aligned_total > 50] = "#808080"
mapping.colors_Corr [anno_mapped_sub$percent_reads_aligned_total <= 50] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("percent_reads_aligned_total > 50", "percent_reads_aligned_total <= 50")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()


NMS <- annoNew[match(anno_mapped_sub$exp_component_name, annoNew$exp_component_name), 
               "marker_sum_norm_label"] 
NMS <- NMS > 0.6
percalign <- anno_mapped_sub$percent_reads_aligned_total > 25
ngenes <- anno_mapped_sub$Genes.Detected > 1000
pass = NMS & percalign & ngenes

jpeg(file.path(mappingFolder,'204_337_AIT115_umap_compoundQC.jpg'), quality = 100)
main="UMAP Visualization of Corr mapped cells"
pad=0.1 # was 0.1
cex=0.4
pch=19
cex.main=1
cex.legend=0.7
xylim <- range(layout)
xlim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-4.5, 0.5)
ylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)

#par(mar=c(0.2,0.7,1.2,0.7), ps=10)
plot(xlim, ylim, type="n", axes=F, frame=F, xlab="UMAP1", ylab="UMAP2")
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="#aaaaaa", lwd=0.25)  

mapping.colors_Corr <- mapping.colors_Corr
mapping.colors_Corr [pass == 'TRUE'] = "#808080"
mapping.colors_Corr [pass == 'FALSE'] = "#FF69B4"

points(layout[,1], layout[,2], col=mapping.colors_Corr, cex=cex, pch=pch)
mtext(side=3, main, cex=cex.main)

labels.u <- c("QC pass", "QC fail")
legend.pos <- "topleft"
legend.text <- as.character(labels.u)
legend.colors <- c("#808080", "#FF69B4")
legend(legend.pos, legend=legend.text, inset=0.03,
       col=legend.colors,
       bty="n", pch=pch, cex=cex.legend)

dev.off()

# Cluster on UMAP
load(paste0(data_dir, "/20230807_RSC-204-337_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20230807_RSC-204-337_macaque_patchseq_star2.7_samp.dat.Rdata"))

query.metadata <- samp.dat
counts      <- cpmR   # Genes are rows, samples are columns

query.counts   <- counts
query.data   <- logCPM(query.counts)

# Put annotations and counts in the same order
query.metadata <- query.metadata[match(colnames(query.data),query.metadata$exp_component_name),] 
rownames(query.metadata) <- query.metadata$exp_component_name  

annoNew <- annoNew[match(rownames(query.metadata),annoNew$exp_component_name),] 
rownames(annoNew) <- annoNew$exp_component_name 

inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
#inds2 = annotations_mapped$library_prep_pass_fail == "Pass"
#inds3 = annotations_mapped$Genes.Detected >= 1000
#inds4 = annotations_mapped$percent_reads_aligned_total >= 50
annoNew = annoNew[inds1,] 

load(file=file.path(mappingFolder,"NHP_BG_204_337_AIT115_umap_roi.Rdata"))
layout <- mapping.umap$layout
install.packages("dbscan")
library("dbscan")
library("scales")
cl <- hdbscan(layout, minPts = 10)
n_cl = length(unique(cl$cluster))
hex_codes <- hue_pal()(n_cl)
hex_codes[4] = hex_codes[1]
hex_codes[1] = "#5A5A5A"
colors = hex_codes[cl$cluster+1]
png(file.path(mappingFolder,'NHP_BG_204_337_AIT115_umap_roi_hdbscan.png'), width = 1000, height = 1000)
plot(layout, col=colors, pch=20, cex = 1.5, xlab = 'UMAP1', ylab = 'UMAP2')
#colors = unique(cl$cluster+1)
legendtxt = c("Cluster #1 (noise)", "Cluster #2", "Cluster #3", "Cluster #4", "Cluster #5", "Cluster #6")
#legendtxt[colors==1] = paste(legendtxt[colors==1], "(noise)", sep = " ")
# first two arguments are x,y positions
legend(-4, 11, legend = legendtxt, 
       fill = hex_codes
)
dev.off()

inds3 = annoNew$Genes.Detected >= 1000
inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
inds6 = annoNew$marker_sum_norm_label >= 0.6
annoNew$compound_qc_pass = inds3 & inds4 & inds6

for (clus in sort(unique(cl$cluster))) {
  if (clus != 0){
    cat("Cluster #:", clus+1, '\n')
    anno_sub = annoNew[cl$cluster == clus,]
    n = dim(anno_sub)[1]
    print("Fraction RNA_amplification_pass")
    print(sum(anno_sub$rna_amplification_pass_fail=="Pass")/n)
    print("Fraction compound QC pass")
    print(sum(anno_sub$compound_qc_pass)/n)
    print("Fraction Tree and Corr subclasses agree")
    print(sum(anno_sub$level3.subclass_Corr == anno_sub$level3.subclass_Tree)/n)
    print("mean NMS")
    print(mean(anno_sub$marker_sum_norm_label))
    type_counts <- table(anno_sub$level3.subclass_Tree)
    print(type_counts)
    species_counts <- table(anno_sub$species)
    print(species_counts)
    acute <- table(anno_sub$cell_specimen_project)
    print(acute)
    contam <- table(anno_sub$contaminationType_label)
    print(contam)
  }
  if (clus == 4){
    clus5_list = anno_sub$cell_name
    write.csv(clus5_list, file.path(mappingFolder,'NHP_BG_204_337_umap_cluster5.csv')) 
  }
}


# Debugging before subsetting
hist_info <- hist(layout[,2]) 
table(mapping.labels_Corr[layout[,2]< -4])
mean(annotations_mapped$score.Corr[layout[,2]< -4])
mean(annotations_mapped$score.Corr[layout[,2]> -4])

mean(annotations_mapped$Genes.Detected[layout[,2]< -4])
mean(annotations_mapped$Genes.Detected[layout[,2]> -4])

hist_info <- hist(annotations_mapped$Genes.Detected) 
hist_info <- hist(annotations_mapped$score.Corr)    #0.2 cutoff?
hist_info <- hist(annotations_mapped$score.Tree).    # start with 0.1?