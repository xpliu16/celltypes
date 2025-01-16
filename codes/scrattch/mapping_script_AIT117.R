# singularity shell --cleanenv docker://jeremyinseattle/scrattch:0.1.1

suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)
library(scrattch.patchseq)
library(scrattch.taxonomy)
library (dplyr)
library(stringr)
#library(scrattch.patchseq)
library(reticulate)
#cell_type_mapper <- import("cell_type_mapper")
#reticulate::use_python("/usr/bin/python3")

source("/home/xiaoping.liu/scrattch/mapping/run_mappings.R") 

run_mappings(refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Macaque/", 
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_117",  
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/",
             data_fn = "202411107_RSC-204-378_macaque_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = "HMBA_Macaque_BG_082024_AIT.h5ad",
             hierarchy <- c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label"),
             #class_colname = 'Neighborhood_label', 
             #neigh_colname = 'Class_label', 
             #subclass_colname = 'Subclass_label',  
             #cluster_colname = 'Group_label',      # HACK because hierarchy is different, to match mouse whole brain
             low_level = 'Group_label',
             proj_strs = "qIVSCC-MET",
             roi_strs = "STR|PALGPi|PALGPe|PAL_GPe|HYSTN|OT_L",
             off_target = c("Immune", "Astro-Epen", "Vascular", "OPC-Oligo"),
             off_target_level = 'Class_label'
)
refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Macaque/" 
mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_117"  
data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/" 
data_fn = "202411107_RSC-204-378_macaque_patchseq_star2.7" 
mode = 'patchseq'                                                                                  
h5ad_fn = "HMBA_Macaque_BG_082024_AIT.h5ad"
hierarchy <- c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label", "Cluster_label")
#class_colname = 'Class_label' 
#neigh_colname = 'Neighborhood_label' 
#subclass_colname = 'Subclass_label'  
low_level = 'Group_label'
#cluster_colname = 'Group_label'      # HACK because hierarchy is different, to match mouse whole brain

proj_strs = "qIVSCC-MET" 
roi_strs = "STR|PALGPi|PALGPe|PAL_GPe|HYSTN|OT_L"
off_target = c("Immune", "Astro-Epen", "Vascular", "OPC-Oligo")
off_target_level = 'Class_label'