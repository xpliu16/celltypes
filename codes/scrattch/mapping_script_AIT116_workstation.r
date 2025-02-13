suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)
library (dplyr)
library(stringr)

source("run_mappings_noHANN.R")

run_mappings(refFolder = "GreatApes_Macaque_NCBI",
             mappingFolder = "mapping/GreatApes_Macaque_NCBI", 
             data_dir =  "R_Object/",
             data_fn = "20250106_RSC-204-380_macaque_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = 'GreatApes_Macaque_NCBI.h5ad', 
             class_colname = 'class_label',
             neigh_colname = 'neighborhood_label',
             subclass_colname = 'subclass_label', 
             cluster_colname = 'cluster_label', 
             proj_strs = "qIVSCC-MET",
             roi_strs = "",
             off_target = "Nonneuron"
)


refFolder = "GreatApes_Macaque_NCBI"
mappingFolder = "mapping/GreatApes_Macaque_NCBI"
data_dir =  "R_Object/"
data_fn = "20250106_RSC-204-380_macaque_patchseq_star2.7"
mode = 'patchseq'
h5ad_fn = 'GreatApes_Macaque_NCBI.h5ad'
class_colname = 'class_label'
neigh_colname = 'neighborhood_label'
subclass_colname = 'subclass_label' 
cluster_colname = 'cluster_label'
proj_strs = "qIVSCC-MET"
roi_strs = ""
off_target = "Nonneuron"

run_mappings(refFolder = "GreatApes_Human",
             mappingFolder = "mapping/GreatApes_Human", 
             data_dir =  "R_Object/",
             data_fn = "20250106_RSC-122-380_human_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = 'GreatApes_Human.h5ad', 
             class_colname = 'class_label',
             neigh_colname = 'neighborhood_label',
             subclass_colname = 'subclass_label', 
             cluster_colname = 'cluster_label', 
             proj_strs = "qIVSCC-MET",
             roi_strs = "",
             off_target = "glia"
)
refFolder = "GreatApes_Human"
mappingFolder = "mapping/GreatApes_Human"
data_dir =  "R_Object/"
data_fn = "20250106_RSC-122-380_human_patchseq_star2.7"
mode = 'patchseq'
h5ad_fn = 'GreatApes_Human.h5ad'
class_colname = 'class_label'
neigh_colname = 'neighborhood_label'
subclass_colname = 'subclass_label'
cluster_colname = 'cluster_label'
proj_strs = "qIVSCC-MET"
roi_strs = ""
off_target = "glia"





# First time to fix non-relative path
AIT.anndata = read_h5ad(file.path(refFolder,h5ad_fn))
AIT.anndata$uns$taxonomyDir=refFolder
write_h5ad(AIT.anndata, file.path(refFolder,h5ad_fn))