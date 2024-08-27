library(scrattch.mapping)

mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116"
refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_116"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

#annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_AIT_116_RSC-204-366_ann_map_full_QC.csv"))
annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_AIT_116_RSC-204-366_sub_QC.csv"))

HANN_obj <- fromJSON(readLines(file.path(mappingFolder, "20240621_RSC-204-366_neurons_results.json")))
HANN_res <- do.call(cbind, HANN_obj$results)
orig_query <- read_h5ad(file.path(mappingFolder, '20240621_RSC-204-366_query.h5ad'))
HANN_cell_names <- orig_query$obs$exp_component_name
HANN_res <- cbind(HANN_cell_names, HANN_res)
rownames(HANN_res) = HANN_cell_names
HANN_res_sub <- HANN_res[annoNew$exp_component_name,]
HANN_res_sub2 <- select(HANN_res_sub, c("level3.subclass.assignment", "cluster.assignment"))
colnames(HANN_res_sub2) <- c("Subclass_HANN", "Cluster_HANN")

patch_anno <- merge(HANN_res_sub2, annoNew, by.x= 'row.names', by.y='exp_component_name') 

proj_strs = "qIVSCC-MET"
roi_strs = "STR"
inds0 = grepl(proj_strs, patch_anno$cell_specimen_project)
inds1 = grepl(roi_strs, patch_anno$roi)
#inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
#inds2 = annotations_mapped$library_prep_pass_fail == "Pass"   # Chucks good samples
inds3 = patch_anno$Genes.Detected >= 1000
inds4 = patch_anno$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
#inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
#inds6 = annotations_mapped$score.Corr > 0.6
inds6 = patch_anno$marker_sum_norm_label >= 0.6

patch_anno_roi_proj = patch_anno[inds0&inds1,]
table(patch_anno_roi_proj$level3.subclass.assignment)

png(file.path(mappingFolder, 'HANN_bootstrapping_probability_hist.png'), width = 500, height = 350)
hist(patch_anno_roi_proj$level3.subclass.bootstrapping_probability)
dev.off()

png(file.path(mappingFolder, 'HANN_avg_correlation_hist.png'), width = 500, height = 350)
hist(patch_anno_roi_proj$level3.subclass.avg_correlation)
dev.off()

patch_anno_sub = patch_anno[inds0&inds1&inds3&inds4&inds6,]
table(patch_anno_sub$level3.subclass.assignment)


# For morphology

annoNew <- read.csv(file.path(mappingFolder, "NHP_BG_AIT_116_RSC-204-366_roi_QC.csv"))

HANN_obj <- fromJSON(readLines(file.path(mappingFolder, "20240621_RSC-204-366_neurons_results.json")))
HANN_res <- do.call(cbind, HANN_obj$results)
orig_query <- read_h5ad(file.path(mappingFolder, '20240621_RSC-204-366_query.h5ad'))
HANN_cell_names <- orig_query$obs$exp_component_name
HANN_res <- cbind(HANN_cell_names, HANN_res)
rownames(HANN_res) = HANN_cell_names
HANN_res_sub <- HANN_res[annoNew$exp_component_name,]
HANN_res_sub2 <- select(HANN_res_sub, c("level3.subclass.assignment", "cluster.assignment"))
colnames(HANN_res_sub2) <- c("Subclass_HANN", "Cluster_HANN")

annoNew <- merge(HANN_res_sub2, annoNew, by.x= 'row.names', by.y='exp_component_name') 

inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
inds3 = annoNew$Genes.Detected >= 1000
inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
inds6 = annoNew$marker_sum_norm_label >= 0.6
annoNew$compound_qc_pass = inds3 & inds4 & inds6
annoNew$BG_ROI = inds1
annoNew$acute = NaN
annoNew$acute[annoNew$cell_specimen_project == "qIVSCC-METa"] = "TRUE"
annoNew$acute[annoNew$cell_specimen_project == "qIVSCC-METc"] = "FALSE"
annoNew$revisit = (annoNew$rna_amplification_pass_fail=="Fail") & (annoNew$compound_qc_pass == TRUE) 

desired_columns = c('cell_name', 'cell_id', 'score.Corr', 'Subclass_Corr', 
                    'score.Tree', 'Subclass_Tree', 'Subclass_HANN', 'rna_amplification_pass_fail', 
                    'compound_qc_pass', 'roi', 'species', 'postPatch_classification', 
                    'acute', 'Virus', 'creCell', 'go_no_go_63x', 'revisit')
# Or striatal ROI?
anno_morpho = annoNew [desired_columns]

write.csv (anno_morpho, file.path(mappingFolder, 'NHP_BG_AIT_116_RSC-204-366_anno_morpho_with_HANN.csv'))
