# Just ran DLPFC mapping
annoNew_DLPFC<- annoNew

# Load up M1 mapping
region = "M1"
refFolder = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_",region)
mode = paste0("CrossAreal_",region,"_patchseq")
h5ad_fn = paste0("CrossAreal_",region,".h5ad")
a <- strsplit(refFolder,'/')[[1]]
taxname <- a[length(a)]
b <- strsplit(data_fn, '_')[[1]]
dataname <- b[2]

vars<-load(file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_full_QC.Rdata', sep='_')))
rownames(annoNew) = annoNew$exp_component_name

# Subset to project
annoNew_DLPFC = annoNew_DLPFC[grepl(proj_strs, annoNew_DLPFC$cell_specimen_project),]
# OR
annoNew_DLPFC <- subset(annoNew_DLPFC, cell_specimen_project %in% c("hIVSCC-MET","hIVSCC-METc"))
annoNew_M1 = annoNew_M1[grepl(proj_strs, annoNew_M1$cell_specimen_project),]

# See waht fraction mapped to same class
dim(annoNew_DLPFC)
dim(annoNew_M1)
sum(annoNew_DLPFC$Class_Corr == annoNew_M1$Class_Corr)
# Surprisingly many cells are class ambiguous, but maybe they are bad quality?

# Subset to good quality by aggressive NMS (considering that the ROI might be wrong)
inds3 = annoNew_M1$Genes.Detected >= 1000
inds4 = annoNew_M1$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
#inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
#inds6 = annotations_mapped$score.Corr > 0.6
inds6 = annoNew_M1$marker_sum_norm_label >= 0.6
annoNew_M1_sub = annoNew_M1[inds3&inds4&inds6,]

inds3 = annoNew_DLPFC$Genes.Detected >= 1000
inds4 = annoNew_DLPFC$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
#inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
#inds6 = annotations_mapped$score.Corr > 0.6
inds6 = annoNew_DLPFC$marker_sum_norm_label >= 0.6
annoNew_DLPFC_sub = annoNew_DLPFC[inds3&inds4&inds6,]

both = intersect(rownames(annoNew_DLPFC_sub), rownames(annoNew_M1_sub))
length(both)
annoNew_M1_inter<-annoNew_M1_sub[both,]
annoNew_DLPFC_inter<-annoNew_DLPFC_sub[both,]
sum(annoNew_DLPFC_inter$Class_Corr == annoNew_M1_inter$Class_Corr)

sum(annoNew_DLPFC_inter$CrossArea_subclass_Corr == annoNew_M1_inter$CrossArea_subclass_Corr)

conf <- annoNew_DLPFC_inter$CrossArea_subclass_Corr != annoNew_M1_inter$CrossArea_subclass_Corr
df <- data.frame(annoNew_DLPFC_inter$CrossArea_subclass_Corr[conf], annoNew_M1_inter$CrossArea_subclass_Corr[conf])
df[1:30,]

conf <- annoNew_DLPFC_inter$CrossArea_subclass_Tree != annoNew_M1_inter$CrossArea_subclass_Tree
df <- data.frame(annoNew_DLPFC_inter$CrossArea_subclass_Tree[conf], annoNew_M1_inter$CrossArea_subclass_Tree[conf])
df[1:30,]


table(annoNew_M1_sub$roi)



#######

data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/"
data_fn = "20240321_RSC-122-359_human_patchseq_star2.7"
load(paste0(data_dir, paste0(data_fn, "_samp.dat.Rdata")))
df <- samp.dat
df2 <- subset(df, cell_specimen_project %in% c("hIVSCC-MET","hIVSCC-METc"))
df2$manualRoi2 <- gsub('[0-6./]', '', df2$manualRoi)
df2$manualRoi2 <- toupper(df2$manualRoi2)
df2$manualRoi2 <- gsub('_L', '', df2$manualRoi2)
unique(df2$manualRoi2)

mappingFolder = "/home/xiaoping.liu/scrattch/mapping/Human_Cortical"
#load(file.path(mappingFolder, "CrossAreal_M1_RSC-122-359_ann_map_full_QC.Rdata"))
load(file.path(mappingFolder, "CrossAreal_MTG_RSC-122-359_ann_map_full_QC.Rdata"))

annoNew <- subset(annoNew, cell_specimen_project %in% c("hIVSCC-MET","hIVSCC-METc"))
table(annoNew$Virus)
table(annoNew[annoNew$Virus=="CN1390",'CrossArea_subclass_Corr'])

#######
# To get passing counts from reload

annoNew <- subset(annoNew, cell_specimen_project %in% c("hIVSCC-MET","hIVSCC-METc"))
#roi_strs = "TCx|Tcx|TEa"
#roi_strs = "OCx|VISp"
roi_strs = "FCx|MOp"

inds1 = grepl(roi_strs, annoNew$roi)
#inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
#inds2 = annotations_mapped$library_prep_pass_fail == "Pass"   # Chucks good samples
inds3 = annoNew$Genes.Detected >= 1000
inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
#inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
#inds6 = annotations_mapped$score.Corr > 0.6
inds6 = annoNew$marker_sum_norm_label >= 0.6

annoNew_sub = annoNew[inds1&inds3&inds4&inds6,]

subclass_colname = "CrossArea_subclass_label" 
cluster_colname = "CrossArea_cluster_label"

type_counts_Corr = table(annoNew_sub[paste0(str_replace(subclass_colname,'_label',''),'_Corr')])
subclass_counts_Tree_M1 = table(annoNew_sub[paste0(str_replace(subclass_colname,'_label',''),'_Tree')])

cluster_counts_Tree_M1 = table(annoNew_sub[paste0(str_replace(cluster_colname,'_label',''),'_Tree')])

save(subclass_counts_Tree_M1, cluster_counts_Tree_M1, subclass_counts_Tree_V1, cluster_counts_Tree_V1, 
     subclass_counts_Tree_MTG, cluster_counts_Tree_MTG, file=file.path(mappingFolder, 'cell_counts.Rdata'))
