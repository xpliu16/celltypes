## Load the relevant libraries
suppressPackageStartupMessages({
  library(DropletUtils)
  library(dplyr)
  library(gplots)
  library(Matrix)
  library(Matrix.utils)   # Not sure if this is needed... I don't think it is
  library(cowplot)
  library(Seurat)
  library(scrattch.hicat) # for logCPM
  library(VENcelltypes)   # For mitochondrial and sex gene lists
  library(future)
  library(ggplot2) 
  library(feather)
  library(mfishtools)
  library(data.table)
  library(Rdpack)
  library(Seurat)
  library(gridExtra)
  library(viridis)
  library(robustbase)     # robustbase outlier
  library(scater)         # scater and robustbase outlier
  library(dendextend)     # For rfTreeMapping
  library(gridExtra)
  library(reticulate)
  library(QCR)
})

anndata = reticulate::import("anndata")
pandas = reticulate::import("pandas")

## My scripts
source("/home/nelson.johansen/scripts/R/taxonomy_helpers.R")
source("/home/nelson.johansen/scripts/R/plotting_helpers.R")
source("/home/nelson.johansen/scripts/R/processing_ATAC.R")
source("/home/nelson.johansen/scripts/R/processing_RNA.R")
source("/home/nelson.johansen/scripts/R/seurat_h5ad.R")

options(stringsAsFactors=FALSE)
options(future.globals.maxSize = 32000 * 1024^2)
options(future.rng.onMisuse="ignore")

data(mito_genes)  
data(sex_genes)  

## -------------
## Read / process all data and metadata
## -----------------------------
setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/NelsonJ/PatchSeq/nhp")

## 
adata = anndata$read_h5ad("./data/nhp_patchseq_exon_intron.h5ad");
data = adata$X; rownames(data) = adata$obs_names$tolist(); colnames(data) = adata$var_names$tolist(); 
meta = as.data.frame(adata$obs);

## Standard Seurat pipeline, will only use logCPM, clustering and UMAP.
rnaseq.data = CreateSeuratObject(t(as.matrix(data)), meta.data=meta)
rnaseq.data = NormalizeData(rnaseq.data, normalization.method = "LogNormalize", scale.factor = 1e6)
rnaseq.data = FindVariableFeatures(rnaseq.data, selection.method = "vst", nfeatures = 2000)
rnaseq.data = ScaleData(rnaseq.data, features = rownames(rnaseq.data))
rnaseq.data = RunPCA(rnaseq.data, features = VariableFeatures(object = rnaseq.data))
rnaseq.data = FindNeighbors(rnaseq.data, dims = 1:30)
rnaseq.data = FindClusters(rnaseq.data, resolution = 0.5)
rnaseq.data = RunUMAP(rnaseq.data, dims = 1:30)

## QC
rnaseq.data = subset(rnaseq.data, cells = colnames(rnaseq.data)[grepl("STR", rnaseq.data$roi)])
rnaseq.data = subset(rnaseq.data, subset=library_prep_pass_fail == "Pass")
rnaseq.data = subset(rnaseq.data, subset=nFeature_RNA >= 1000)
rnaseq.data = subset(rnaseq.data, subset=percent_reads_aligned_total >= 50)

## -------------
## Perform labeling using basic methods (VERY OLD)
## -----------------------------

## Store results
mappingOut = list()

## Load in a reference for which to to basic cell type labeing (more advanced labeling done seperatly)
GEXRef = loadGEXRef(hGenes = rownames(rnaseq.data), refFolder='/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104', nGenes=2000)

## This section will perform initial mapping using three methods: Correlation-based, tree-based, and Seurat based, and will calculate some QC metrics.
mappingOut[["AIT11_3"]] = taxonomy_mapping(rnaseq.data, GEXRef, tree.map=FALSE, seurat.map=FALSE, dims=30, k.weight=15, label.cols = c("cluster_label","class_label","subclass_label", "supertype_label"))

##
rnaseq.data@meta.data = cbind(rnaseq.data@meta.data, mappingOut[["AIT11_3"]])

## -------------
## Filter using meta-cell QC approaches
## -----------------------------
rnaseq.data@meta.data$confident_annotation = rnaseq.data@meta.data$score.Corr >= 0.5
rnaseq.data@meta.data$meta_cell_conf_anno = QC_cluster_class(rnaseq.data, "confident_annotation", class.cutoff=0.8)
rnaseq.data@meta.data$meta_cell_class = QC_cluster_class(rnaseq.data, "class_Corr", class.cutoff=0.8)

## 
rnaseq.data = subset(rnaseq.data, subset=meta_cell_conf_anno == TRUE)
rnaseq.data = subset(rnaseq.data, subset=meta_cell_class == TRUE)
rnaseq.data = subset(rnaseq.data, subset=subclass_Corr != "NN")

## Filter out clusters with poor meta-cell QC / mapping 
rnaseq.data = subset(rnaseq.data, cells=colnames(rnaseq.data)[!rnaseq.data@meta.data$seurat_clusters %in% c(1,3,4,5)])

## Write mapping and QC info table to csv
meta.data = rnaseq.data@meta.data[,c("patched_cell_container", "class_Corr", "subclass_Corr", "supertype_Corr", "score.Corr", "structure", "roi", "postPatch_classification")]
write.table(meta.data, file="nhp_patchseq_mapping_AIT112.csv", row.names=F, sep=",")

## Save .h5ad to get back into python
seurat2h5ad(rnaseq.data, file.name="nhp_patchseq_processed_exon_intron")

## -------------
## Simple visualization
## -----------------------------
for(anno in c("class_Corr", "subclass_Corr", "supertype_Corr")){
  order = names(sort(table(meta.data[,anno])))
  plot.me = meta.data[,anno,drop=F]
  plot.me[,anno] = factor(plot.me[,anno], order)
  anno_bar = ggplot(plot.me, aes_string(x=anno, fill=anno)) + geom_bar() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(anno_bar, file=paste0(anno, "_anno.png"), dpi=600)
}

## Figures
features = c("DRD1", "DRD2", "RXFP1", "PVALB", "SST", "TAC3", "VIP", "CHAT")

pdf("./figures/rdg.pdf", width=14, height=14)
RidgePlot(rnaseq.data, features = features, group.by = 'subclass_Corr', ncol = 3)
dev.off()

pdf("./figures/vln.pdf")
VlnPlot(rnaseq.data, features = features, group.by = 'subclass_Corr', ncol = 1)
dev.off()

pdf("./figures/features.pdf", width=14, height=14)
FeaturePlot(rnaseq.data, features = features)
dev.off()

pdf("./figures/umap.pdf")
DimPlot(rnaseq.data, group.by='supertype_Corr')
dev.off()
