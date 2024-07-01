library(anndata)
library("scrattch.mapping")

mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116"
data_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

data_batch = "20240520_RSC-204-363"
## Load in data
# cpmR
var <- load(file.path(data_path, paste0(data_batch, '_macaque_patchseq_star2.7_cpm.Rdata')))
#samp.dat
var <- load(file.path(data_path, paste0(data_batch, '_macaque_patchseq_star2.7_samp.dat.Rdata')))
rownames(samp.dat) <- samp.dat$cell_id    # or exp_component_name or cell_name

## Add mat and samp.dat to new .h5ad and save
AIT115_anndata <- AnnData(
  X = Matrix::t(logCPM(cpmR)),    # Matrix of observations x features
  obs = samp.dat     # Metadata for observations
)
##
write_h5ad(AIT115_anndata, file.path(mappingFolder, paste0(data_batch, '_query.h5ad')))