#refFolder <- "/home/xiaoping.liu/scrattch/reference/NHP_BG_AIT_114"
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115"    # MAKE THIS MORE SPECIFIC
panelFolder <- "/home/xiaoping.liu/scrattch/MERFISH_panel"
dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)
library (dplyr)

# This remakes taxonomy from feather files and reloads, very slow, use read_h5ad instead for the time being
AIT.anndata <- loadTaxonomy(refFolder)

#AIT.anndata <- read_h5ad(file.path(refFolder,"AIT_114_taxonomy.h5ad"))
#annoReference   = feather(file.path(refFolder,"anno.feather")) 
#exprReference   = feather(file.path(refFolder,"data.feather"))
#annoReference = as.data.frame(annoReference[match(exprReference[[sample_id]], annoReference[[sample_id]]),])
#rownames(annoReference) = rownames(datReference) = annoReference[[sample_id]]

#AIT.anndata$uns$dend = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114/reference.rda"

## Add in the off.target annotation.
AIT.anndata$obs$off_target = AIT.anndata$obs$level1.class_label

## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "level3.subclass_label", ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                    class.column = "off_target", ## The column by which off-target types are determined.
                                    off.target.types = c("NN"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = refFolder)

AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")

# Was "/20220907_RSC-204-310_macaque_patchseq_star2.7_cpm.Rdata"
"/20221215_RSC-204-318_macaque_patchseq_star2.7_cpm.Rdata"
load(paste0(data_dir, "/20231219_RSC-204-350_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20231219_RSC-204-350_macaque_patchseq_star2.7_samp.dat.Rdata"))

#  dim(samp.dat) - 1551  153.  (data frame)
#  names(samp.dat)
#  dim(cpmR) - 17133  1551

query.metadata <- samp.dat
counts      <- cpmR   # Genes are rows, samples are columns

query.counts   <- counts
query.data   <- logCPM(query.counts)

# Put annotations and counts in the same order
query.metadata <- query.metadata[match(colnames(query.data),query.metadata$exp_component_name),] 
rownames(query.metadata) <- query.metadata$exp_component_name  

## Get marker genes from dendrogram
#load(AIT.anndata$uns$dend)
#dend = reference$dend
dend <- readRDS(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
#dend_file = paste0(refFolder,'/dend.RData')
#dend <- readRDS(dend_file)
allMarkers = unique(unlist(get_dend_markers(dend)))

## Alternative get genes from MERFISH panel
#panel_df = read.csv(file.path(panelFolder, "AIT_115_MERFISH_panel.xlsx"))


## Subset query data to just those markers
query.data = query.data[intersect(rownames(query.data), allMarkers),]

# Confirm
#query.metadata$exp_component_name[1:5]
#colnames(counts)[1:5]

# Filter out unclassifieds
#kp          <- annotations$broad_type!="Unclassified"
#counts      <- counts[,kp]
#annotations <- annotations[kp,]

## Identify the offtarget cell types manually.
print(unique(AIT.anndata$obs$level1.class_label))

#AIT.anndata$uns$dend$patchseq = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115/patchseq/dend.RData"
# For testing small batch:
#query.datasub = query.data[,1:10]

query.mapping <- taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.data, 
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = FALSE, 
#                                  label.cols = c("class_label", "subclass_label", "supertype_label", "cluster_label")  # Currently just use subclass_label and cluster. Columns to map against
                                  label.cols = c("level1.class_label", "level2.neighborhood_label", "level3.subclass_label", "cluster_label")
)

#label.cols = c("class", "subclass", "supertype_label", "cluster") 

# HACK because I didn't map against supertype (in AIT114 this is subclass)
#query.mapping$supertype_Corr <- AIT.anndata$obs$supertype_label[match(query.mapping$cluster_Corr, AIT.anndata$obs$cluster_label)]
#query.mapping$supertype_Tree <- AIT.anndata$obs$supertype_label[match(query.mapping$cluster_Tree, AIT.anndata$obs$cluster_label)]
#query.mapping$level3_subclass_Corr <- AIT.anndata$obs$level3.subclass_label[match(query.mapping$cluster_Corr, AIT.anndata$obs$cluster_label)]
#query.mapping$level3_subclass_Tree <- AIT.anndata$obs$level3.subclass_label[match(query.mapping$cluster_Tree, AIT.anndata$obs$cluster_label)]

write.csv(query.mapping, file.path(mappingFolder,"NHP_BG_204_350_AIT115_mapping_withglia.csv"), row.names=FALSE)
save(query.mapping, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_mapping_withglia.Rdata"))

#write.csv(query.mapping, file.path(mappingFolder,"NHP_BG_204_350_AIT115_mapping.csv"), row.names=FALSE)
#save(query.mapping, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_mapping.Rdata"))

# Variable renaming
#clusters  <- unique(query.mapping$cluster)   
clusters <- unique(AIT.anndata$uns$clusterInfo$cluster_label)
#subclass_levels <- unique(query.mapping$subclass_Corr)
#subclass_levels <- unique(AIT.anndata$uns$clusterInfo$subclass_label)

annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = 0, all=TRUE) 
#annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = "exp_component_name", all=TRUE) 
# Alt hack:
#rownames(query.mapping) = rownames(query.metadata)
#annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = 0, all=TRUE) 
# annotations_mapped$clusters <- factor(annotations_mapped$cluster_Tree, levels=clusters)
# annotations_mapped$class   <- annotations_mapped$class_Tree
# annotations_mapped$class[!is.element(annotations_mapped$class,c("IN","MSN"))] = "Non-neuronal"
annotations_mapped <- annotations_mapped[match(rownames(query.metadata),annotations_mapped$exp_component_name),]   # merge resorts things
rownames(annotations_mapped) <- annotations_mapped$exp_component_name
#type_counts_Corr = table(annotations_mapped$supertype_Corr)
#type_counts_Tree = table(annotations_mapped$supertype_Tree)
type_counts_Corr = table(annotations_mapped$level3.subclass_Corr)
type_counts_Tree = table(annotations_mapped$level3.subclass_Tree)

write.csv(annotations_mapped, file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full_withglia.csv"), row.names=FALSE)
save(annotations_mapped, type_counts_Corr, type_counts_Tree, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full_withglia.Rdata"))

#dir.create(file.path(mappingFolder, "NHP_BG_RSC_204_350_map_full"))

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = file.path(mappingFolder, "NHP_BG_RSC_204_350_map_full"),
                      query.data     = counts,  # Don't need log-normalized data here
                      query.metadata = query.metadata,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = TRUE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

dir(file.path(mappingFolder, "NHP_BG_RSC_204_350_map_full"))

# leaf node names, each also has an ID, 0 appears to be unclassified, 
# In this taxonomy there are 49 types - max(annotations$primary_type_id)

# Don't need to run this if you are running buildMappingDirectory, but will do it to generate Rdata and csv version
#writePatchseqQCmarkers(counts = query.counts,
#                       metadata = query.metadata,
#                       subsample = 100,  # Default of 100 is reasonable
#                       subclass.column = "level3.subclass_Corr",  # default
#                       class.column = "level1.class_Corr",  # default
#                       #off.target.types = "Non-neuron",  # default is various iterations of non-neuronal
#                       off.target.types = c("NN"),
#                       num.markers = 50,     # Default of 50 is probably fine
#                       shinyFolder = paste0(refFolder,"/")    # don't need to paste0 if your refFolder path has /
#)

# Don't need to run this if you are running buildMappingDirectory, but will do it to generate Rdata and csv version
annoNew = applyPatchseqQC (AIT.anndata = AIT.anndata, ## A patchseq taxonomy object.
                           query.data = query.counts, ## Counts are required here.
                           #query.data = query.counts,
                           query.metadata = annotations_mapped, ## Results of the previous mapping or AIT.anndata$obs, no mapping is required.
                           #query.metadata = annotations_mapped,
                           verbose=FALSE)
# Ran the contents of the function directly in R to avoid error in /R/patchseq_output.R
save(annoNew, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full_QC.Rdata"))
write.csv(annoNew, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full_QC.csv"))

# annotations_mapped$cluster <- factor(annotations_mapped$clusters, levels=clusters)  # Make into discrete levels
inds1 = grepl("STR",annotations_mapped$roi)
#inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
#inds2 = annotations_mapped$library_prep_pass_fail == "Pass"   # Chucks good samples
inds3 = annoNew$Genes.Detected >= 1000
inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
#inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
#inds6 = annotations_mapped$score.Corr > 0.6
inds6 = annoNew$marker_sum_norm_label >= 0.6

#query.counts   <- counts
#query.data   <- logCPM(query.counts)

#query.counts_sub <- counts[,inds1&inds3&inds4&inds6]
#query.metadata_sub <- query.metadata[inds1&inds3&inds4&inds6,]
#query.mapping_sub <- query.mapping[inds1&inds3&inds4&inds6,]

#query.counts_sub_df <- as.data.frame(query.counts_sub)
# Don't need if you ran buildPatchseqTaxonomy
#writePatchseqQCmarkers(counts = query.counts_sub,
#                       metadata = query.metadata_sub,
#                       subsample = 100,  # Default of 100 is reasonable
#                       subclass.column = "level3.subclass_Corr",  # default
#                       class.column = "level1.class_Corr",  # default
#                       #off.target.types = "Non-neuron",  # default is various iterations of non-neuronal
#                       off.target.types = c("NN"),
#                       num.markers = 50,     # Default of 50 is probably fine
#                       shinyFolder = paste0(refFolder,"/")    # don't need to paste0 if your refFolder path has /
#)
# counts = paste0(refFolder,"/counts.feather"),  # can also read matrix directly
# metadata = paste0(refFolder,"/anno.feather"),  # can also read matrix directly
# AIT.anndata$uns$QC_markers<-paste0(refFolder,"/QC_markers.RData")
# NHP_QC <- load(paste0(refFolder,"/patchseq/QC_markers.RData")) # Just to take a look... looks reasonable

# Don't need to run this is you are running buildMappingDirectory
#annoNew = applyPatchseqQC (AIT.anndata = AIT.anndata, ## A patchseq taxonomy object.
#                           query.data = query.counts_sub, ## Counts are required here.
#                           #query.data = query.counts,
#                           query.metadata = anno_mapped_sub, ## Results of the previous mapping or AIT.anndata$obs, no mapping is required.
#                           #query.metadata = annotations_mapped,
#                           verbose=FALSE)
# Ran the contents of the function directly in R to avoid error in /R/patchseq_output.R

annoNew_roi = annoNew[inds1,]
save(annoNew_roi, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC.Rdata"))
write.csv(annoNew_roi, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC.csv"))

# To merge in mappings to full taxonomy (including glia)
annoNew_roi_with_glia = annotations_mapped[inds1,]
var<-load(file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC.Rdata"))
annoNew_roi_with_glia = annoNew_roi_with_glia[,2:11]
colnames(annoNew_roi_with_glia) <- paste0(colnames(annoNew_roi_with_glia), '_wglia')
annotations_mapped2 <- merge(x = annoNew_sub, y = annoNew_roi_with_glia, by.x = "exp_component_name", by.y = 0, all=TRUE) 
annoNew_wglia <- annotations_mapped2[match(rownames(annoNew_roi_with_glia),annotations_mapped2$exp_component_name),]   # merge resorts things
save(annoNew_wglia, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC_withglia.Rdata"))
write.csv(annoNew_wglia, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC_withglia.csv"))

annoNew_sub = annoNew[inds1&inds3&inds4&inds6,]
save(annoNew_sub, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_sub_QC.Rdata"))
write.csv(annoNew_sub, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_sub_QC.csv"))

# Sampling counts
load(file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_sub_QC.Rdata"))

dim(annoNew_sub)
type_counts_Corr = table(annoNew_sub$level3.subclass_Corr)
type_counts_Tree = table(annoNew_sub$level3.subclass_Tree)

main_subclasses = c('D1-Matrix', 'D2-Matrix', 'D1-Striosome', 'D2-Striosome', 'D1-ShellOT', 'D2-ShellOT', 'D1D2-Hybrid', 'D2-Hybrid-MCHR2', 'D1-NUDAP', 'PVALB-COL19A1-ST18', 'SST_Chodl', 'CCK-FBXL7', 'CHAT', 'CCK-VIP-TAC3', 'LHX6-TAC3-PLPP4', 'TAC3-LHX8-PLPP4')
remaining_subclasses = setdiff(names(type_counts_Tree),main_subclasses)
main_subclasses <- c(main_subclasses, remaining_subclasses)
type_counts_Tree <- data.frame(type_counts_Tree)
#type_counts_Tree$Var1 <- factor(type_counts_Tree$Var1, levels = main_subclasses)
type_counts_Tree$Var1 <- ordered(type_counts_Tree$Var1, levels = main_subclasses)

colorkey <- read.csv(file = file.path(mappingFolder, 'colortable.csv'))
#colors = list()
rownames(type_counts_Tree) = type_counts_Tree$Var1
type_counts_Tree['color'] <- NA
for (c in type_counts_Tree$Var1) {
  tmp = colorkey['colz'][colorkey['subclass']==c]
  if (length(tmp)==0) {
    tmp = '#000000'
  } # Do we need this?
  #colors <- append (colors, tmp)
  type_counts_Tree[c,'color'] <- tmp
  print(type_counts_Tree[c,'color'])
  #print(colorkey$colz[colorkey['subclass']==c])
}

png(file.path(mappingFolder,'NHP_BG_AIT115_sampling_counts.png'), width = 500, height = 350)
tmp <- par("mar")
tmp[1] = tmp[1]+7
par(mar = tmp)
#barplot(names = type_counts_Tree$Var1, height = type_counts_Tree$Freq, las=2, col=type_counts_Tree$color)
ggplot(type_counts_Tree,aes(x= Var1, y = Freq)) +
       geom_bar(stat= 'Identity', fill = type_counts_Tree$color) +
       #scale_fill_manual(values=type_counts_Tree$color) +
       xlab("") +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_blank()) +
       geom_hline(yintercept=10,linetype=2, color = 'gray')
dev.off()

# Striatal subclasses only (at least 5% of all cells are in dSTR or vSTR)
subclass = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
             "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
             "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
             "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
             "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
             "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")  

# cumsum
df_NMS <- subset(annoNew_roi, level3.subclass_Tree %in% c('SN_STH', 'SLC17A7-SATB2', 'MEIS2', 'D1-Matrix', 'D2-Matrix', 'PVALB-COL19A1-ST18'))
png(file.path(mappingFolder,'NHP_BG_AIT115_cumNMS_SN_STH.png'), width = 500, height = 250)
ggplot(df_NMS, aes(marker_sum_norm_label, color=level3.subclass_Tree)) + 
  stat_ecdf(geom = "step", size=1) + 
  labs(x="NMS", y="Cumulative Fraction") +
  scale_y_continuous(breaks=seq(0,1,0.1), labels = seq(0,100,10)) + 
  #scale_color_manual(values=c("#96ceb4", "#ff6f69", "#ffcc5c", "#90697c"), name='Subclass') + 
  scale_color_manual(values=c("#758f0b", "#ff6f69", "#ffcc5c", "#90697c", "#007c8f", "#cf0690"), name='Subclass') +
  theme_minimal()
dev.off()

postpatch1 = annoNew_sub['postPatch_classification'][annoNew_sub['level3.subclass_Tree'] == 'SN_STH']
table(postpatch1)/length(postpatch1)
postpatch2 = annoNew_sub['postPatch_classification'][annoNew_sub['level3.subclass_Tree'] == 'D1-Matrix']
table(postpatch2)/length(postpatch2)

# Count cells with ephys
df_ephys = read.csv(file=file.path("NHP_ephys_features_20231207.csv"))
df_id = read.csv("custom_report_20231204.csv")

df2 = merge(annoNew_sub, df_id, by.x='cell_name', by.y='cell_specimen_name.', all.x = FALSE, all.y = FALSE)
# checked these numbers against python by merging against annoNew_roi instead
df3 = merge(df2, df_ephys, by.x ='cell_specimen_id.', right_on='cell_name', all.x = FALSE, all.y = FALSE)
sum(is.na(df3['tau']))
sum(is.na(df3['adapt_hero']))
sum(is.na(df3['upstroke_downstroke_ratio_short_square']))
grep("sag", colnames(df3))
df_ephys = df3[265:358]
df_ephys <- df_ephys[,colSums(is.na(df_ephys))<nrow(df_ephys)]  # drop any rows with all Na; nothing dropped

df_short = subset(type_counts_Tree, select = -c(color))
df_short['source'] = 'T'
type_counts_ephys = table(df3$level3.subclass_Tree)
type_counts_ephys <- data.frame(type_counts_ephys)
type_counts_ephys['source'] = 'E'
rownames(type_counts_ephys) = type_counts_ephys$Var1
type_counts = rbind(df_short,type_counts_ephys)
type_counts$source <- ordered(type_counts$source, levels = c('T','E'))

png(file.path(mappingFolder,'NHP_BG_AIT115_sampling_counts_wephys.png'), width = 650, height = 350)
tmp <- par("mar")
tmp[1] = tmp[1]+7
par(mar = tmp)
#barplot(names = type_counts_Tree$Var1, height = type_counts_Tree$Freq, las=2, col=type_counts_Tree$color)
ggplot(type_counts,aes(x= Var1, y = Freq, fill = source)) +
  geom_bar(stat= 'Identity',position="dodge") +
  #scale_fill_manual(values=type_counts_Tree$color) +
  xlab("") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=10,linetype=2, color = 'gray') +
  scale_fill_manual(values = c("#62b7c4", "#d455ac"))
dev.off()


dir.create(file.path(mappingFolder, "NHP_BG_204_350_AIT115_map_sub_patchseqQC"))
#unlink(file.path(mappingFolder, "NHP_BG_204_350_AIT115_map_sub_patchseqQC/dend.RData"))

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = file.path(mappingFolder, "NHP_BG_204_350_AIT115_map_sub_patchseqQC"),
                      query.data     = query.counts_sub,  # Don't need log-normalized data here
                      query.metadata = query.metadata_sub,
                      query.mapping  = query.mapping_sub,
                      doPatchseqQC   = TRUE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

unique(annotations_mapped$roi)

anno_mapped_roi = annoNew[inds1,]
query.data_roi = query.data[,inds1]
write.csv(anno_mapped_roi, file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi.csv"), row.names=FALSE)
type_counts_Corr_roi = table(anno_mapped_roi$level3.subclass_Corr)
type_counts_Tree_roi = table(anno_mapped_roi$level3.subclass_Tree)

#anno_mapped_sub = annotations_mapped[inds1&inds3&inds4&inds6,]
query.data_sub = query.data[,inds1&inds3&inds4&inds6]

rownames(annoNew) <- annoNew$exp_component_name   # After running apply_PatchseqQC
anno_mapped_sub <- annoNew[inds1&inds3&inds4&inds6,]  

type_counts_Corr_QC = table(anno_mapped_sub$level3.subclass_Corr)
type_counts_Tree_QC = table(anno_mapped_sub$level3.subclass_Tree)
write.csv(anno_mapped_sub, file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_sub_QC.csv"), row.names=FALSE)

type_counts_Corr_QC_df = as.data.frame(type_counts_Corr_QC, row.names = NULL,
              responseName = "Freq")
# Optionally remove any non-neuronal types
subclasses = unique(AIT.anndata$obs$level3.subclass_label[AIT.anndata$obs$level1.class_label!="NN"])
type_counts_Corr_QC_df <- filter(type_counts_Corr_QC_df, Var1 %in% subclasses)

# Sort by predominant region, print out up to top 3 locations
roi_df <- data.frame(matrix(ncol = 11, nrow = 0))

for (type in subclasses) {
  t <- table(AIT.anndata$obs$roi_label[AIT.anndata$obs$level3.subclass_label == type])
  roi_df <- rbind(roi_df, t)
  colnames(roi_df) <- dimnames(t)[[1]]
  rownames(roi_df)[length(rownames(roi_df))] <- type
}
roi_df <- roi_df/rowSums(roi_df)
roi_df_grouped <- data.frame(matrix(ncol = 0, nrow = nrow(roi_df)))
dSTR = c("Macaque CaB", "Macaque CaH", "Macaque CaT", "Macaque PuC", "Macaque PuPV",  
        "Macaque PuR")
roi_df_grouped$dSTR = rowSums(roi_df[, dSTR])
vSTR = "Macaque NAC"
roi_df_grouped$vSTR = roi_df[, vSTR]
GPe = "Macaque GPe"   
roi_df_grouped$GPe = roi_df[, GPe]
GPi = "Macaque GPi"
roi_df_grouped$GPi = roi_df[, GPi]
SN_VTA = "Macaque SN-VTA"
roi_df_grouped$SN_VTA = roi_df[, SN_VTA]
STH = "Macaque STH"
roi_df_grouped$STH = roi_df[, STH]
rownames(roi_df_grouped) <- rownames(roi_df)
roi_df_sorted <- roi_df_grouped[order(-roi_df_grouped$dSTR, -roi_df_grouped$vSTR), ]

jpeg(file.path(mappingFolder,'204_350_AIT115_roi_distr_all.jpg'), quality = 100, width = 1000, height = 2000)
plot.new()
par(mfrow = c(nrow(roi_df_sorted), 1))  # Divide the plotting area into multiple rows
par(oma = c(2, 30, 3, 1))
# Generate bar graphs for each row
for (i in 1:nrow(roi_df_sorted)) {
  row <- roi_df_sorted[i, ]
#  if (i == 1){
#    main = 
#  } else {
#    main = NULL
#  }
  if (i == nrow(roi_df_sorted)) {
    xaxt = "s"
  } else {
    xaxt = "n"
  }
  par(mar = c(1.8, 1, 1.8, 10))
  par(mgp = c(3, 1.6, 0))    # Offset the x-axis labels
  barplot(t(as.matrix(row)), names.arg = rownames(row), horiz = TRUE, 
          las = 1, cex.axis = 2.8, cex.names = 2.8,
          xaxt = xaxt, 
          col=c("red","purple","green","darkgreen","cyan","black")
          )
  if (i ==1) {
    par(mar = c(1.8, 1, 3.4, 10))
    title(main = "Fraction of cells in each region", line = 1.8, cex.main =2.5)
    legend(x = 1.02, y = 1.9, legend=colnames(roi_df_sorted), cex=1.3, xpd=TRUE,
           fill = c("red","purple","green","darkgreen","cyan","black"), 
           box.lwd = 0, y.intersp=0.8)
  }
}

dev.off()

roi_df_sorted$STR = roi_df_sorted$dSTR + roi_df_sorted$vSTR
str_types = rownames(roi_df_sorted[roi_df_sorted$STR > 0.05,])

str_types = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
              "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
              "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
              "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
              "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
              "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
              "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   


type_counts_Tree_QC_df = as.data.frame(type_counts_Tree_QC, row.names = NULL,
                                       responseName = "Freq")
# Optionally remove any non-neuronal types - should no longer be necessary
# n_subclasses = unique(AIT.anndata$obs$level3.subclass_label[AIT.anndata$obs$level1.class_label!="NN"])
type_counts_Tree_QC_df <- filter(type_counts_Tree_QC_df, Var1 %in% n_subclasses)

print(type_counts_Corr_QC_df, row.names = FALSE)
print(type_counts_Tree_QC_df, row.names = FALSE)

write.csv(annotations_mapped, file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_QC.csv"), row.names=FALSE)
save(annotations_mapped, anno_mapped_sub, type_counts_Corr, type_counts_Tree, type_counts_Corr_QC,
     type_counts_Tree_QC, type_counts_Corr_roi, type_counts_Tree_roi, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115ann_map_QC.Rdata"))

D1_expr <- query.data_sub["DRD1",]
D1_expr_D1_types = D1_expr[is.element(anno_mapped_sub$subclass_Tree,c("D1-Matrix","D1-ShellOT","D1-Striosome", "D1-ICj"))]
D1_expr_D1_types <- na.omit(D1_expr_D1_types)
length(D1_expr_D1_types)
D1_expr_D1_types = mean(D1_expr_D1_types)
D1_expr_D2_types = D1_expr[is.element(anno_mapped_sub$subclass_Tree,c("D2-Matrix","D2-ShellOT","D2-Matrix / D2-Striosome", "D2-Striosome"))]
D1_expr_D2_types <- na.omit(D1_expr_D2_types)
D1_expr_D2_types = mean(D1_expr_D2_types)

D2_expr <- query.data_sub["DRD2",]
D2_expr_D1_types = D2_expr[is.element(anno_mapped_sub$level3.subclass_Tree,c("D1-Matrix","D1-ShellOT","D1-Striosome", "D1-ICj"))]
D2_expr_D1_types <- na.omit(D2_expr_D1_types)
length(D2_expr_D1_types)
D2_expr_D1_types = mean(D2_expr_D1_types)
D2_expr_D2_types = D2_expr[is.element(anno_mapped_sub$level3.subclass_Tree,c("D2-Matrix","D2-ShellOT","D2-Matrix / D2-Striosome", "D2-Striosome"))]
D2_expr_D2_types <- na.omit(D2_expr_D2_types)
D2_expr_D2_types = mean(D2_expr_D2_types)

RXFP1_expr <- query.data_sub["RXFP1",]
RXFP1_expr_D1D2_types = RXFP1_expr[is.element(anno_mapped_sub$level3.subclass_Tree,c("D1D2-Hybrid"))]
#RXFP1_expr_D1D2_types <- na.omit(RXFP1_expr_D1D2_types)
length(RXFP1_expr_D1D2_types)
RXFP1_expr_D1D2_types = mean(RXFP1_expr_D1D2_types)
# May exclude "D2-Matrix / D2-Striosome" which may express RXFP1, but doesn't make a difference
indsA = is.element(anno_mapped_sub$level3.subclass_Tree, c("D1-Matrix","D1-ShellOT","D1-Striosome", "D2-Matrix / D2-Striosome", "D1-ICj", "D2-Matrix","D2-ShellOT", "D2-Striosome"))
indsB = !is.element(anno_mapped_sub$cluster_Tree, c("53_MSN", "54_MSN", "55_MSN"))  # THESE MAY NOT BE UP TO DATE WITH AIT 11.5
RXFP1_expr_other_types = RXFP1_expr[indsA&indsB]
#RXFP1_expr_other_types <- na.omit(RXFP1_expr_other_types)
length(RXFP1_expr_other_types)
RXFP1_expr_other_types = mean(RXFP1_expr_other_types)

# Filter by confidence
inds5 <- anno_mapped_sub$score.Tree >= 0.6
query.data_sub2 = query.data_sub[,inds5]
anno_mapped_sub2 = anno_mapped_sub[inds5,]
RXFP1_expr <- query.data_sub2["RXFP1",]
RXFP1_expr_D1D2_types = RXFP1_expr[is.element(anno_mapped_sub2$subclass_Tree,c("D1D2 Hybrid"))]
#RXFP1_expr_D1D2_types <- na.omit(RXFP1_expr_D1D2_types)
length(RXFP1_expr_D1D2_types)
RXFP1_expr_D1D2_types = mean(RXFP1_expr_D1D2_types)
# May exclude "D2-Matrix / D2-Striosome" which may express RXFP1, but doesn't make a difference
indsA = is.element(anno_mapped_sub2$subclass_Tree, c("D1-Matrix","D1-ShellOT","D1-Striosome", "D2-Matrix / D2-Striosome", "D1-ICj", "D2-Matrix","D2-ShellOT", "D2-Striosome"))
indsB = !is.element(anno_mapped_sub2$cluster_Tree, c("53_MSN", "54_MSN", "55_MSN"))
RXFP1_expr_other_types = RXFP1_expr[indsA&indsB]
#RXFP1_expr_other_types <- na.omit(RXFP1_expr_other_types)
length(RXFP1_expr_other_types)
RXFP1_expr_other_types = mean(RXFP1_expr_other_types)

inds5 <- anno_mapped_sub$score.Corr >= 0.6
query.data_sub2 = query.data_sub[,inds5]
anno_mapped_sub2 = anno_mapped_sub[inds5,]
RXFP1_expr <- query.data_sub2["RXFP1",]
RXFP1_expr_D1D2_types = RXFP1_expr[is.element(anno_mapped_sub2$subclass_Corr,c("D1D2 Hybrid"))]
#RXFP1_expr_D1D2_types <- na.omit(RXFP1_expr_D1D2_types)
length(RXFP1_expr_D1D2_types)
RXFP1_expr_D1D2_types = mean(RXFP1_expr_D1D2_types)
# May exclude "D2-Matrix / D2-Striosome" which may express RXFP1, but doesn't make a difference
indsA = is.element(anno_mapped_sub2$subclass_Corr, c("D1-Matrix","D1-ShellOT","D1-Striosome", "D2-Matrix / D2-Striosome", "D1-ICj", "D2-Matrix","D2-ShellOT", "D2-Striosome"))
indsB = !is.element(anno_mapped_sub2$cluster_Corr, c("53_MSN", "54_MSN", "55_MSN"))
RXFP1_expr_other_types = RXFP1_expr[indsA&indsB]
#RXFP1_expr_other_types <- na.omit(RXFP1_expr_other_types)
length(RXFP1_expr_other_types)
RXFP1_expr_other_types = mean(RXFP1_expr_other_types)

PVALB_expr <- query.data_sub["PVALB",]
df_pv <- data.frame(anno_mapped_sub$subclass_Tree, PVALB_expr)
means <- aggregate(PVALB_expr ~ anno_mapped_sub.subclass_Tree, data = df_pv, FUN = mean)
inds_pv = AIT.anndata$uns$clusterInfo$subclass_label=="IN_str-PVALB"
PVALB_expr_tax = AIT.anndata$X[,colnames(AIT.anndata$X)=="PVALB"]
pv_expr_pv_type_tax = mean(PVALB_expr_tax[inds_pv])
num_above_mean <- sum(PVALB_expr > pv_expr_pv_type_tax, na.rm=TRUE)
num_below_mean <- sum(PVALB_expr < pv_expr_pv_type_tax, na.rm=TRUE)
sub = anno_mapped_sub[PVALB_expr > pv_expr_pv_type_tax,] 
table(sub$subclass_Tree)


D1_expr <- query.data_sub["DRD1",]
D1_expr_D1_types = D1_expr[is.element(anno_mapped_sub$subclass_Corr,c("D1-Matrix","D1-ShellOT","D1-Striosome", "D1-ICj"))]
D1_expr_D1_types <- na.omit(D1_expr_D1_types)
length(D1_expr_D1_types)
D1_expr_D1_types = mean(D1_expr_D1_types)
D1_expr_D2_types = D1_expr[is.element(anno_mapped_sub$subclass_Corr,c("D2-Matrix","D2-ShellOT","D2-Matrix / D2-Striosome", "D2-Striosome"))]
D1_expr_D2_types <- na.omit(D1_expr_D2_types)
D1_expr_D2_types = mean(D1_expr_D2_types)

D2_expr <- query.data_sub["DRD2",]
D2_expr_D1_types = D2_expr[is.element(anno_mapped_sub$subclass_Corr,c("D1-Matrix","D1-ShellOT","D1-Striosome", "D1-ICj"))]
D2_expr_D1_types <- na.omit(D2_expr_D1_types)
length(D2_expr_D1_types)
D2_expr_D1_types = mean(D2_expr_D1_types)
D2_expr_D2_types = D2_expr[is.element(anno_mapped_sub$subclass_Corr,c("D2-Matrix","D2-ShellOT","D2-Matrix / D2-Striosome", "D2-Striosome"))]
D2_expr_D2_types <- na.omit(D2_expr_D2_types)
D2_expr_D2_types = mean(D2_expr_D2_types)

RXFP1_expr <- query.data_sub["RXFP1",]
RXFP1_expr_D1D2_types = RXFP1_expr[is.element(anno_mapped_sub$subclass_Corr,c("D1D2 Hybrid"))]
#RXFP1_expr_D1D2_types <- na.omit(RXFP1_expr_D1D2_types)
length(RXFP1_expr_D1D2_types)
RXFP1_expr_D1D2_types = mean(RXFP1_expr_D1D2_types)
# May exclude "D2-Matrix / D2-Striosome" which may express RXFP1, but doesn't make a difference
indsA = is.element(anno_mapped_sub$subclass_Corr, c("D1-Matrix","D1-ShellOT","D1-Striosome", "D2-Matrix / D2-Striosome", "D1-ICj", "D2-Matrix","D2-ShellOT", "D2-Striosome"))
indsB = !is.element(anno_mapped_sub$cluster_Corr, c("53_MSN", "54_MSN", "55_MSN"))
RXFP1_expr_other_types = RXFP1_expr[indsA&indsB]
#RXFP1_expr_other_types <- na.omit(RXFP1_expr_other_types)
length(RXFP1_expr_other_types)
RXFP1_expr_other_types = mean(RXFP1_expr_other_types)

KCNC1_expr <- query.data_sub["KCNC1",]
KCNC1_expr_COL19A1_types = KCNC1_expr[is.element(anno_mapped_sub$subclass_Corr,c("IN_str-LHX6-COL19A1"))]
length(KCNC1_expr_COL19A1_types)
KCNC1_expr_COL19A1_types = mean(KCNC1_expr_COL19A1_types)
KCNC1_expr_MSN_types = KCNC1_expr[is.element(anno_mapped_sub$class_Corr, c("MSN"))]
length(KCNC1_expr_MSN_types)
KCNC1_expr_MSN_types = mean(KCNC1_expr_MSN_types)
indsA = is.element(anno_mapped_sub$class_Corr, c("IN"))
indsB = !is.element(anno_mapped_sub$subclass_Corr, c("IN_str-LHX6-COL19A1"))
KCNC1_expr_otherIN_types = KCNC1_expr[indsA&indsB]
length(KCNC1_expr_otherIN_types)
KCNC1_expr_otherIN_types = mean(KCNC1_expr_otherIN_types)

KCNC2_expr <- query.data_sub["KCNC2",]
KCNC2_expr_COL19A1_types = KCNC2_expr[is.element(anno_mapped_sub$subclass_Corr,c("IN_str-LHX6-COL19A1"))]
length(KCNC2_expr_COL19A1_types)
KCNC2_expr_COL19A1_types = mean(KCNC2_expr_COL19A1_types)
KCNC2_expr_MSN_types = KCNC2_expr[is.element(anno_mapped_sub$class_Corr, c("MSN"))]
length(KCNC2_expr_MSN_types)
KCNC2_expr_MSN_types = mean(KCNC2_expr_MSN_types)
indsA = is.element(anno_mapped_sub$class_Corr, c("IN"))
indsB = !is.element(anno_mapped_sub$subclass_Corr, c("IN_str-LHX6-COL19A1"))
KCNC2_expr_otherIN_types = KCNC2_expr[indsA&indsB]
length(KCNC2_expr_otherIN_types)
KCNC2_expr_otherIN_types = mean(KCNC2_expr_otherIN_types)

KCNC3_expr <- query.data_sub["KCNC3",]
KCNC3_expr_COL19A1_types = KCNC3_expr[is.element(anno_mapped_sub$subclass_Corr,c("IN_str-LHX6-COL19A1"))]
length(KCNC3_expr_COL19A1_types)
KCNC3_expr_COL19A1_types = mean(KCNC3_expr_COL19A1_types)
KCNC3_expr_MSN_types = KCNC3_expr[is.element(anno_mapped_sub$class_Corr, c("MSN"))]
length(KCNC3_expr_MSN_types)
KCNC3_expr_MSN_types = mean(KCNC3_expr_MSN_types)
indsA = is.element(anno_mapped_sub$class_Corr, c("IN"))
indsB = !is.element(anno_mapped_sub$subclass_Corr, c("IN_str-LHX6-COL19A1"))
KCNC3_expr_otherIN_types = KCNC3_expr[indsA&indsB]
length(KCNC3_expr_otherIN_types)
KCNC3_expr_otherIN_types = mean(KCNC3_expr_otherIN_types)

table(anno_mapped_sub$Level3.subclass_Tree)

#Hack
#rownames(query.mapping) = rownames(query.metadata)

#counts = counts[intersect(rownames(counts), allMarkers),]

# QC files for Rachel
# Optional load annoNew from another run
# Run on server:
load(file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full_QC.Rdata"))
# Run on laptop:
#load(file="/Users/xiaoping.liu/celltypes/NHP_BG_anal/NHP_BG_AIT_115/204_350/NHP_BG_204_350_AIT115_ann_map_full_QC.Rdata")
inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
inds3 = annoNew$Genes.Detected >= 1000
inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
inds6 = annoNew$marker_sum_norm_label >= 0.6
annoNew$compound_qc_pass = inds3 & inds4 & inds6
annoNew$BG_ROI = inds1
annoNew$acute = NaN
annoNew$acute[annoNew$cell_specimen_project == "qIVSCC-METa"] = "TRUE"
annoNew$acute[annoNew$cell_specimen_project == "qIVSCC-METc"] = "FALSE"
annoNew$revisit1 = (annoNew$rna_amplification_pass_fail=="Fail") & (annoNew$compound_qc_pass == TRUE) 

desired_columns = c('exp_component_name', 'cell_name', 'cell_id', 'level3.subclass_Corr', 
                    'level3.subclass_Tree', 'rna_amplification_pass_fail', 'compound_qc_pass', 
                    'BG_ROI', 'roi', 'species', 'postPatch_classification', 'acute', 'Virus', 
                    'creCell')
# Or striatal ROI?
anno_morpho = annoNew[desired_columns]
write.csv(anno_morpho, file.path(mappingFolder,"NHP_BG_204_350_AIT115_anno_morpho.csv"))

annoNew_hiconf = annoNew[annoNew$score.Corr >= 0.6,]
mean(annoNew_hiconf$quality_score_label)
annoNew_loconf = annoNew[annoNew$score.Corr < 0.6,]
mean(annoNew_loconf$quality_score_label)

#df <- read_feather(file.path(mappingFolder, "NHP_BG_204_350_AIT115_map_sub_patchseq_roi/anno.feather"))
load(file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC.Rdata"))
df <- as.data.frame(annoNew_sub)
df_hiconf = df[df$score.Corr_label >= 0.6,]
mean(df_hiconf$quality_score_label)
df_loconf = df[df$score.Corr_label < 0.6,]
mean(df_loconf$quality_score_label)

df_IN <- df[df$class_Corr_label=="IN",]
mean(df_IN$quality_score_label)
df_MSN <- df[df$class_Corr_label=="MSN",]
mean(df_MSN$quality_score_label)
df_NN <- df[df$class_Corr_label=="NN",]
mean(df_NN$quality_score_label)

mean(df_IN$contam_sum_label)
mean(df_MSN$contam_sum_label)
mean(df_NN$contam_sum_label)

load(file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_umap_roi.Rdata"))
layout <- mapping.umap$layout
df_sort <- df[match(rownames(mapping.umap$layout),df$exp_component_name_label),] 
library(RColorBrewer)
#pal <- colorRamp(c("blue", "red"))
nicereds <- brewer.pal(5, "Reds")
reds <- colorRamp(nicereds)
mapping.colors <- reds(df_sort$quality_score_label)

layout1 = layout[,1]
layout2 = layout[,2]
quality = df_sort$quality_score_label
mapping_umap <- data.frame(layout1, layout2, quality)
mapping_umap$contam_sum = df_sort$contam_sum_label
mapping_umap$NMS = df_sort$Norm_Marker_Sum.0.4_label

library("ggplot2")

jpeg(file.path(mappingFolder,'204_350_NMS_hist.jpg'), quality = 100)
hist(df_sort$marker_sum_norm_label) 
dev.off()
df_sort$Norm_Marker_Sum.0.6_calc<-df_sort$marker_sum_norm_label>0.6

jpeg(file.path(mappingFolder,'204_350_umap_quality_score.jpg'), quality = 100)

ggplot(mapping_umap,aes(x=layout1,y=layout2,col=quality))+geom_point()
dev.off()

jpeg(file.path(mappingFolder,'204_350_umap_contam_sum.jpg'), quality = 100)
ggplot(mapping_umap,aes(x=layout1,y=layout2,col=contam_sum))+geom_point()
dev.off()

jpeg(file.path(mappingFolder,'204_350_umap_NMS_Pass_Fail.jpg'), quality = 100)
ggplot(mapping_umap,aes(x=layout1,y=layout2,col=NMS))+geom_point()
dev.off()

# Looking at result of failing neurons with NMS threshold 0.6 instead of 0.4
mapping_umap$NMS = df_sort$Norm_Marker_Sum.0.6_calc
jpeg(file.path(mappingFolder,'204_350_umap_NMS_06.jpg'), quality = 100)
ggplot(mapping_umap,aes(x=layout1,y=layout2,col=NMS))+geom_point()
dev.off()

# Looking at distributions of NMS for nucleated, partial, and no seal
df_nuc = df_sort[df_sort$postPatch_classification_label == "Nucleated",]
df_pn = df_sort[df_sort$postPatch_classification_label == "Partial-Nucleus",]
df_ns = df_sort[df_sort$postPatch_classification_label == "No-Seal",]

jpeg(file.path(mappingFolder,'204_350_NMS_dist_by_patch.jpg'), quality = 100)
par(mfrow=c(3,1))
bins = seq(from = 0, to = max(df_sort$marker_sum_norm_label)+0.1, by = 0.1)
hist(df_nuc$marker_sum_norm_label, 
      main = 'NMS values "Nucleated" samples',
      freq = TRUE, plot = TRUE, breaks=bins)
hist(df_pn$marker_sum_norm_label, 
     main = 'NMS values "Partial-Nucleus" samples',
     freq = TRUE, plot = TRUE, breaks=bins)
hist(df_ns$marker_sum_norm_label, 
     main = 'NMS values "No-Seal" samples',
     freq = TRUE, plot = TRUE, breaks=bins)
dev.off()

jpeg(file.path(mappingFolder,'204_350_quality_hist.jpg'), quality = 100)
hist(quality) 
dev.off()

cor(df_sort$quality_score_label, df_sort$score.Corr_label, method = c("spearman"))
cor(df_sort$quality_score_label, df_sort$score.Tree_label, method = c("spearman"))
cor(df_sort$marker_sum_norm_label, df_sort$score.Corr_label, method = c("spearman"))

# Looking at distributions of NMS for nucleated, partial, and no seal
df_nuc = quality[df_sort$postPatch_classification_label == "Nucleated"]
df_pn = quality[df_sort$postPatch_classification_label == "Partial-Nucleus"]
df_ns = quality[df_sort$postPatch_classification_label == "No-Seal"]

jpeg(file.path(mappingFolder,'204_350_quality_dist_by_patch.jpg'), quality = 100)
par(mfrow=c(3,1))
bins = seq(from = 0, to = max(quality)+0.1, by = 0.1)
hist(df_nuc, 
     main = 'Quality values "Nucleated" samples',
     freq = TRUE, plot = TRUE, breaks=bins)
hist(df_pn, 
     main = 'Quality values "Partial-Nucleus" samples',
     freq = TRUE, plot = TRUE, breaks=bins)
hist(df_ns, 
     main = 'Quality values "No-Seal" samples',
     freq = TRUE, plot = TRUE, breaks=bins)
dev.off()

jpeg(file.path(mappingFolder,'204_350_ngenes_vs_quality.jpg'), quality = 100)
ggplot(df_sort,aes(x=Genes.Detected_label,y=quality_score_label,col=postPatch_classification_label))+geom_point()
dev.off()

jpeg(file.path(mappingFolder,'204_350_ngenesCPM_vs_quality.jpg'), quality = 100)
ggplot(df_sort,aes(x=Genes.Detected.CPM_label,y=quality_score_label,col=postPatch_classification_label))+geom_point()
dev.off()

jpeg(file.path(mappingFolder,'204_350_ngenes_vs_quality_microglia.jpg'), quality = 100)
ggplot(df_sort,aes(x=Genes.Detected_label,y=quality_score_label,col=subclass_Corr_label))+geom_point()
dev.off()

df<-anno_mapped_sub
#viruses = unique(df$Virus_label)
viruses = unique(df$Virus)
viruses = viruses[!is.element(viruses,c("","ZZ_Missing", "None", NA))]
virus_desc = c(
  "CN3738" = "unknown",
  "CN2421" = "DRD2",
  "CN1839" = "unknown",
  "CN1390" = "dlx2.0 pan gabaergic",
  "CN3445" = "SLC5A7(choline transporter)"
)
for (v in viruses) {
  df2 = df[df$Virus%in%v & df$creCell%in%"Positive",]
  mapped_subclass_v = df2$level3.subclass_Corr
  print(paste0(v,': '))
  print(mapped_subclass_v)
}

MALAT1_expr <- query.data_sub["MALAT1",]
MALAT1_expr_Nuc_types = MALAT1_expr[is.element(df_sorted$postPatch_classification_label,c("Nucleated"))]
length(MALAT1_expr_Nuc_types)
MALAT1_expr_Nuc_types = mean(MALAT1_expr_Nuc_types)
MALAT1_expr_MSN_types = MALAT1_expr[is.element(anno_mapped_sub$class_Corr, c("MSN"))]
length(MALAT1_expr_MSN_types)
MALAT1_expr_MSN_types = mean(MALAT1_expr_MSN_types)
indsA = is.element(anno_mapped_sub$class_Corr, c("IN"))
indsB = !is.element(anno_mapped_sub$subclass_Corr, c("IN_str-LHX6-COL19A1"))
MALAT1_expr_otherIN_types = MALAT1_expr[indsA&indsB]
length(MALAT1_expr_otherIN_types)
MALAT1_expr_otherIN_types = mean(MALAT1_expr_otherIN_types)

NEAT1_expr <- query.data_sub["NEAT1",]
MALAT1_expr_Nuc_types = MALAT1_expr[is.element(df_sorted$postPatch_classification_label,c("Nucleated"))]
length(MALAT1_expr_Nuc_types)
MALAT1_expr_Nuc_types = mean(MALAT1_expr_Nuc_types)
MALAT1_expr_MSN_types = MALAT1_expr[is.element(anno_mapped_sub$class_Corr, c("MSN"))]
length(MALAT1_expr_MSN_types)
MALAT1_expr_MSN_types = mean(MALAT1_expr_MSN_types)
indsA = is.element(anno_mapped_sub$class_Corr, c("IN"))
indsB = !is.element(anno_mapped_sub$subclass_Corr, c("IN_str-LHX6-COL19A1"))
MALAT1_expr_otherIN_types = MALAT1_expr[indsA&indsB]
length(MALAT1_expr_otherIN_types)
MALAT1_expr_otherIN_types = mean(MALAT1_expr_otherIN_types)


# Do tabulation of quality by subclass using aggregation

dir.create(file.path(mappingFolder, "NHP_BG_RSC_204_350_map_full"))

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = file.path(mappingFolder, "NHP_BG_RSC_204_350_map_full"),
                      query.data     = counts,  # Don't need log-normalized data here
                      query.metadata = query.metadata,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = TRUE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

dir(file.path(mappingFolder, "NHP_BG_RSC_204_350_map_full"))


anno_hybrid_tree = anno_mapped_sub[anno_mapped_sub$subclass_Tree=='D1D2 Hybrid',]


library(feather)
df <- read_feather(file.path(mappingFolder,"anno.feather"))

# Mean of each gene across all samples
logCPM_mean <- apply(query.data, 1, mean)
# Sort ascending
logCPM_mean2 = sort(-logCPM_mean)
logCPM_mean2[1:5]
#  RBFOX1       DLG2      AGGF1      GAPDH      NRXN3 
#-10.833781 -10.305143  -9.949140  -9.863655  -9.826239 

df_pass = anno_mapped_sub[anno_mapped_sub$library_prep_pass_fail == "Pass",]
df_fail = anno_mapped_sub[anno_mapped_sub$library_prep_pass_fail == "Fail",]

## Analyze prevalences

AIT.anndata <- read_h5ad(file.path(refFolder,"NHP_BG_AIT115_complete.h5ad"))
data = AIT.anndata$X
subclass = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
             "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
             "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
             "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
             "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
             "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   
#subclass_tx = AIT.anndata$obs$level3.subclass_label
#subclass_tx = AIT.anndata$obs$Subclass   # These are for the _complete taxonomy (but has since changed?)
subclass_tx = AIT.anndata$obs$level3.subclass
#region_tx = AIT.anndata$obs$roi_label
#region_tx = AIT.anndata$obs$Region
region_tx = AIT.anndata$obs$roi
region_str = c("Macaque CaB", "Macaque CaH", "Macaque CaT", "Macaque PuC", "Macaque PuPV",  
                      "Macaque PuR", "Macaque NAC")

n = {}
n_norm = {}
inds2 <-is.element(region_tx,region_str)
for (sc in subclass){
  inds <-is.element(subclass_tx,sc)
  n[sc] = dim(data[inds&inds2,])[1]
  n_norm[sc] = n[sc]/n["D1-Matrix"]
  print(sc)
  print(n_norm[sc])
}
save(n_norm, file=file.path(mappingFolder,"NHP_BG_AIT115_complete_striatal_n_norm.Rdata"))
#load(file.path(mappingFolder,"NHP_BG_AIT115_complete_striatal_n_norm.Rdata"))

n_D1 = 172
n_exp = {}
for (sc in subclass) {
  n_exp[sc] = round(n_norm[sc] * n_D1)
}

# Make pie chart of STR Inh types in dorsal striatum
str_types = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
              "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
              "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
              "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
              "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
              "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
              "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   
subclass_inh = c("LHX6-TAC3-PLPP4", "PVALB-COL19A1-ST18", "LHX6-SATB1", 
                 "CCK-VIP-TAC3", "CCK-FBXL7", "SST-RSPO2", "SST_Chodl", "CHAT", 
                 "TAC3-LHX8-PLPP4", "MEIS2", "LHX6-LHX8-GBX1", "LHX6_SST", 
                 "SST-ADARB2", "SLC17A6") 
dSTR = c("Macaque CaB", "Macaque CaH", "Macaque CaT", "Macaque PuC", "Macaque PuPV",  
         "Macaque PuR")
n = {}
n_norm_inh = {}
inds2 <-is.element(region_tx,dSTR)
for (sc in subclass_inh){
  print(sc)
  print(n[sc])
  inds <-is.element(subclass_tx,sc)
  n[sc] = dim(data[inds&inds2,])[1]
  #n_norm[sc] = n[sc]/n["D1-Matrix"]
}
save(n, file=file.path(mappingFolder,"NHP_BG_AIT115_complete_striatal_n_inh.Rdata"))
n_norm_inh = n/sum(n)

# Pie Chart with Percentages
pct <- round(n_norm_inh * 100, digits = 1)
jpeg(file.path(mappingFolder,'NHP_BG_AIT115_dSTR_inh_proportions.jpg'), quality = 100)
lbls <- names(n_norm_inh)
lbls <- paste(lbls, pct)
# add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(n,labels = lbls, col=rainbow(length(lbls)), radius = 1, cex = 0.3,
    main="Inhibitory subclass proportions (dorsal striatum) AIT11.5")
dev.off()

install.packages("ggforce")
library(ggforce)
#library(ggrepel)
rpie = 1
rlabel = 0.6 
n2 <- sort(n)
n2 <- n2[!(names(n2) %in% c('SST-ADARB2', 'LHX6_SST'))]
n2 <- c(n2[11:12], n2[1:10])
dat <- data.frame(n2)
#dat <- dat[!(names(n2) %in% c('SST-ADARB2', 'LHX6_SST')),] # 3 and 4 cells
n_norm_inh = n2/sum(n2)
pct <- round(n_norm_inh * 100, digits = 1)
lbls <- names(n_norm_inh)
lbls <- paste(lbls, pct)
lbls <- paste(lbls,"%",sep="")
dat$Group = names(n2)

dat$end_angle = 2*pi*cumsum(pct/100)
dat$start_angle = lag(dat$end_angle, default = 0)
dat$mid_angle = 0.5*(dat$start_angle + dat$end_angle) 

rlabel = 1.05 * rpie # now we place labels outside of the pies

dat <- mutate(dat,
              hjust = ifelse(mid_angle>pi, 1, 0),
              vjust = ifelse(mid_angle<pi/2 | mid_angle>3*pi/2, 0, 1))

png(file.path(mappingFolder,'NHP_BG_AIT115_dSTR_inh_proportions2.png'), width = 2000, height = 2000)

pie.4 <-
  ggplot(dat) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = rainbow(length(lbls)))) +
  geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = lbls,
                hjust = hjust, vjust = vjust), size = 5) +
  coord_fixed() +
  scale_x_continuous(limits = c(-2.2, 2.2), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-2.2, 2.2), name = "", breaks = NULL, labels = NULL) +
  #facet_grid(.~Group, switch = "both")+
  #theme_void()+
  #scale_fill_grey()+
  theme(legend.position='none') +
  #theme(legend.position="bottom", legend.direction="vertical", legend.margin = margin(30, 0, 10, 0))+
  #theme(plot.title = element_text(size = 12, hjust = 0.5))+
  labs(title = "Inhibitory Subclass Proportions (dorsal striatum) AIT 11.5") +
  guides (fill =  guide_legend (title.theme = element_text (face = "bold")))

pie.4

dev.off()


D1-Matrix          D2-Striosome             D2-Matrix 
1.000000000           0.184523917           1.170428466 
D2-Hybrid-MCHR2           D1D2-Hybrid       LHX6-TAC3-PLPP4 
0.062019987           0.252218227           0.134954913 
D1-Striosome         SLC17A7-SATB2    PVALB-COL19A1-ST18 
0.237309576           0.021775665           0.173464409 
LHX6-SATB1          CCK-VIP-TAC3             CCK-FBXL7 
0.004300920           0.021125106           0.045286156 
SST-RSPO2             SST_Chodl            D2-ShellOT 
0.007969351           0.064911361           0.235755462 
CHAT            D1-ShellOT              D1-NUDAP 
0.038093861           0.234743481           0.162639825 
TAC3-LHX8-PLPP4                 MEIS2 SN_STH_GPe-MEIS2-OTX2 
0.016245911           0.098198312           0.003975640 
LHX6-LHX8-GBX1              LHX6_SST        NAc-CCK-SEMA3A 
0.051032763           0.000993910           0.029220955 
GP-LHX6            SST-ADARB2               SLC17A6 
0.008439200           0.012920831           0.013264181 
WDR49-ADAM12                D1-ICj              NAc-LHX8 
0.006595948           0.060285162           0.022697291 


D1-Matrix          D2-Striosome             D2-Matrix 
197                    36                   231 
D2-Hybrid-MCHR2           D1D2-Hybrid       LHX6-TAC3-PLPP4 
12                    50                    27 
D1-Striosome         SLC17A7-SATB2    PVALB-COL19A1-ST18 
47                     4                    34 
LHX6-SATB1          CCK-VIP-TAC3             CCK-FBXL7 
1                     4                     9 
SST-RSPO2             SST_Chodl            D2-ShellOT 
2                    13                    46 
CHAT            D1-ShellOT              D1-NUDAP 
8                    46                    32 
TAC3-LHX8-PLPP4                 MEIS2 SN_STH_GPe-MEIS2-OTX2 
3                    19                     1 
LHX6-LHX8-GBX1              LHX6_SST        NAc-CCK-SEMA3A 
10                     0                     6 
GP-LHX6            SST-ADARB2               SLC17A6 
2                     3                     3 
WDR49-ADAM12                D1-ICj              NAc-LHX8 
1                    12                     4 
