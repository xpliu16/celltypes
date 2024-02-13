refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
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

## Add in the off.target annotation.
AIT.anndata$obs$off_target = AIT.anndata$obs$level1.class_label

## Setup the taxonomy for patchseqQC to infer off.target contamination
#AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
#                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
#                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
#                                    subclass.column = "level3.subclass_label", ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
#                                    class.column = "off_target", ## The column by which off-target types are determined.
#                                    off.target.types = c("NN"), ## The off-target class.column labels for patchseqQC.
#                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
#                                    taxonomyDir = refFolder)

AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")

load(paste0(data_dir, "/20231219_RSC-204-350_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20231219_RSC-204-350_macaque_patchseq_star2.7_samp.dat.Rdata"))

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

## Identify the offtarget cell types manually.
print(unique(AIT.anndata$obs$level1.class_label))

query.mapping <- taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.data, 
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = FALSE, 
                                  label.cols = c("level1.class_label", "level2.neighborhood_label", "level3.subclass_label", "cluster_label")
)

write.csv(query.mapping, file.path(mappingFolder,"NHP_BG_204_350_AIT115_mapping.csv"), row.names=FALSE)
save(query.mapping, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_mapping.Rdata"))

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

write.csv(annotations_mapped, file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full.csv"), row.names=FALSE)
save(annotations_mapped, type_counts_Corr, type_counts_Tree, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_full.Rdata"))

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
save(annoNew_wglia, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC.Rdata"))
write.csv(annoNew_wglia, file=file.path(mappingFolder,"NHP_BG_204_350_AIT115_ann_map_roi_QC.csv"))

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