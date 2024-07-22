## Load scrattch.mapping
library(scrattch.mapping)

## Load in example count data
library(tasic2016data)

library(Seurat)

setwd('/home/xiaoping.liu/scrattch')
bins = 0:38 / 2

#jpeg('tasic_2016_log_counts_dist.jpg')
jpeg('log_counts_dists.jpg')
par(mfrow=c(4,1)) 

tasic_counts <- tasic_2016_counts
query.logCPM <- logCPM(tasic_counts)
hist_info <- hist(query.logCPM[query.logCPM!=0], 
                  main = "Tasic 2016 log2CPM (dropping zeroes)",
                  freq = TRUE, plot = TRUE, breaks=bins)
#dev.off()

# Compare BG 11.4 taxonomy
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114"
AIT.anndata <- read_h5ad(file.path(refFolder,"AIT_114_taxonomy.h5ad"))
ann_counts <- AIT.anndata$X
#jpeg('204-324_ann_log_counts_dist.jpg')
ann_counts_sub <- ann_counts[1:10000,]
hist_info <- hist(ann_counts_sub[ann_counts_sub != 0], 
                  main = "NHP_BG_taxonomy X values (dropping zeroes)", 
                  freq = TRUE, plot = TRUE, breaks=bins)
#dev.off()

# BG raw data
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"
load(paste0(data_dir, "/20230309_RSC-204-324_macaque_patchseq_star2.7_cpm.Rdata"))
counts <- cpmR 
log_counts <- logCPM(counts)
#jpeg('204-324_log_counts_dist.jpg')
hist_info <- hist(log_counts[log_counts!=0], 
                  main = "BG Patch-seq data after log2CPM (dropping zeroes)",
                  freq = TRUE, plot = TRUE, breaks=bins)
log2(max(counts))
max(log_counts)
# Looks like logCPM does normalize, but if you put in CPM, it isn't affected by unneeded normalization
#log_counts <- logCPM(counts*4)
#max(log_counts)
#dev.off()

# BG raw data, processed a la Nelson
load(paste0(data_dir, "/20230309_RSC-204-324_macaque_patchseq_star2.7_samp.dat.Rdata"))
meta <- samp.dat
rnaseq.data = CreateSeuratObject(counts, meta.data=meta)   # Not sure if you want to transform, but doing it to get cells on the rows
rnaseq.data = NormalizeData(rnaseq.data, normalization.method = "LogNormalize", scale.factor = 1e6)
rnaseq.data = FindVariableFeatures(rnaseq.data, selection.method = "vst", nfeatures = 2000)
rnaseq.data = ScaleData(rnaseq.data, features = rownames(rnaseq.data))
#cts = GetAssayData(rnaseq.data, slot = "counts")
#log_cts1 = GetAssayData(rnaseq.data, slot = "data")
#OR
log_cts1 = as.matrix(rnaseq.data@assays$RNA@data)
#log_cts1 = as.matrix(rnaseq.data@assays$RNA@scale.data)
# max(log_cts1 from Seurat matches computed log_cts1 below)
max(counts)
counts_norm <- sweep(counts,2,colSums(counts),`/`)
# After dividing columns by sum and multiplying by a million, you have CPM
log_cts2 <- log(counts_norm*1e6)
max(log_cts2)

hist_info <- hist(log_cts1[log_cts1!=0], 
                  main = "BG Patch-seq data after log e CPM (dropping zeroes)",
                  freq = TRUE, plot = TRUE, breaks=bins)


dev.off()

