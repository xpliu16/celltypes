#refFolder <- "/home/xiaoping.liu/scrattch/reference/NHP_BG_AIT_114"
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping"
#dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

library(umap)
library(scrattch.mapping)

load(file.path(mappingFolder,"NHP_BG_RSC_204_324_mapping.Rdata"))
load(file.path(mappingFolder,"NHP_BG_RSC_204_324_ann_map.Rdata"))
#load(file.path(refFolder,"NHP_BG_RSC_204_324_mapping.Rdata"))

load(paste0(data_dir, "/20230309_RSC-204-324_macaque_patchseq_star2.7_samp.dat.Rdata"))
load(paste0(data_dir, "/20230309_RSC-204-324_macaque_patchseq_star2.7_cpm.Rdata"))

query.metadata <- samp.dat
counts      <- cpmR   # Genes are rows, samples are columns

query.counts   <- counts
query.data   <- logCPM(query.counts)

# Put annotations and counts in the same order
query.metadata <- query.metadata[match(colnames(query.data),query.metadata$exp_component_name),] 
rownames(query.metadata) <- query.metadata$exp_component_name  

#Hack
rownames(query.mapping) = rownames(query.metadata)

AIT.anndata <- read_h5ad(file.path(refFolder,"AIT_114_taxonomy.h5ad"))
# Assign levels
clusters <- unique(AIT.anndata$uns$clusterInfo$cluster)
subclass_lavels <- unique(AIT.anndata$uns$clusterInfo$subclass_label)

annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = 0, all=TRUE) 
annotations_mapped <- annotations_mapped[match(rownames(query.metadata),annotations_mapped$exp_component_name),] 
rownames(annotations_mapped) <- annotations_mapped$exp_component_name 

#annotations_mapped$cluster <- factor(query.mapping$cluster, levels=clusters)  # Make into discrete levels
inds1 = grepl("STR",annotations_mapped$roi)
inds2 = annotations_mapped$library_prep_pass_fail == "Pass"
inds3 = annotations_mapped$Genes.Detected >= 1000
inds4 = annotations_mapped$percent_reads_aligned_total >= 50
anno_mapped_sub = annotations_mapped[inds1&inds2&inds3&inds4,]
query.data_sub = query.data[,inds1&inds2&inds3&inds4]
# check n cells

#query.counts_t <- t(query.counts)
#mapping.data <- query.counts_t
query.data_t <- t(query.data_sub)
mapping.data <- query.data_t
mapping.labels_Corr <- anno_mapped_sub$subclass_Corr
mapping.labels_Tree <- anno_mapped_sub$subclass_Tree
# check rows are aligned - yes

#mapping.umap <- umap(mapping.data)
#save(mapping.umap, file=file.path(mappingFolder,"NHP_BG_RSC_204_324_umap_sub_scoreCorr.Rdata"))
load(file=file.path(mappingFolder,"NHP_BG_RSC_204_324_umap_sub.Rdata"))

# Optional: Further subset by score.Corr
inds5 <- anno_mapped_sub$score.Corr >= 0.6
mapping.umap$layout <- mapping.umap$layout[inds5,]
anno_mapped_sub <- anno_mapped_sub[inds5,]
mapping.labels_Corr <- mapping.labels_Corr[inds5]
mapping.labels_Tree <- mapping.labels_Tree[inds5]

mapping.colors_Tree <- 
  AIT.anndata$uns$clusterInfo$subclass_color[match(mapping.labels_Tree, AIT.anndata$uns$clusterInfo$subclass_label)]

mapping.colors_Corr <- 
  AIT.anndata$uns$clusterInfo$subclass_color[match(mapping.labels_Corr, AIT.anndata$uns$clusterInfo$subclass_label)]

layout <- mapping.umap$layout

jpeg(file.path(mappingFolder,'204_324_umap_Tree.jpg'), quality = 100)
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