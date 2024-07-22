library(Seurat)
library(patchwork)
library(scrattch.mapping)

refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/seurat_integration/NHP_BG_AIT_115"    # MAKE THIS MORE SPECIFIC
dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

adata = loadTaxonomy(refFolder)
data = adata$X; rownames(data) = adata$obs_names$tolist(); colnames(data) = adata$var_names$tolist() 
meta = as.data.frame(adata$obs)

