# singularity shell --cleanenv docker://jeremyinseattle/scrattch:0.7.1

source("/home/xiaoping.liu/scrattch/mapping/run_mappings_AIT117.R") 

run_mappings(refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Macaque/", 
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_117",  
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/",
             data_fn = "20250106_RSC-204-380_macaque_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = "HMBA_Macaque_BG_112024_AIT_v2.h5ad",
             hierarchy <- c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label"),
             class_colname = 'Neighborhood_label', 
             #neigh_colname = 'Class_label', 
             subclass_colname = 'Subclass_label',  
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
data_fn = "20250106_RSC-204-380_macaque_patchseq_star2.7" 
mode = 'patchseq'                                                                                  
#h5ad_fn = "HMBA_Macaque_BG_082024_AIT.h5ad"
h5ad_fn = "HMBA_Macaque_BG_112024_AIT_v2.h5ad"
hierarchy <- c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label", "Cluster_label")
class_colname = 'Class_label' 
#neigh_colname = 'Neighborhood_label' 
subclass_colname = 'Subclass_label'  
low_level = 'Group_label'
#cluster_colname = 'Group_label'      # HACK because hierarchy is different, to match mouse whole brain
cluster_colname = 'Cluster_label'
proj_strs = "qIVSCC-MET" 
roi_strs = "STR|PALGPi|PALGPe|PAL_GPe|HYSTN|OT_L"
off_target = c("Immune", "Astro-Epen", "Vascular", "OPC-Oligo")
off_target_level = 'Class_label'

mode = 'patchseq-str' 
off_target = c("Immune", "Astro-Epen", "Vascular", "OPC-Oligo", "F M Glut", "MB Dopaminergic", 
"MB-GABA", "BG GABA Glut", "GP PVALB GABA", "SN PVALB GABA", "TH MEIS2 GABA", "SN STH GABA", 
"BN MEIS2 GABA", "CTX-CGE GABA", "CTX-MGE GABA", "NDB SI LHX8 GABA", "GP STRv LHX8 GABA", 
"GP SOX6 KCNA1 GABA", "NDB SI LHX6 LHX8 GBX1")



a <- strsplit(refFolder,'/')[[1]]
taxname <- a[length(a)]
b <- strsplit(data_fn, '_')[[1]]
dataname <- b[2]
vars<-load(file=file.path(mappingFolder,paste(taxname, dataname, 'sub_QC.Rdata', sep='_')))
vars<-load(file=file.path(mappingFolder, paste(taxname, dataname, 'roi_QC.Rdata', sep='_')))
vars<-load(file=file.path(mappingFolder, paste(taxname, dataname, 'roi_proj_QC.Rdata', sep='_')))
vars<-load(file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_full_QC.Rdata', sep='_')))

# To merge in mappings to full taxonomy (including glia)
annoNew_roi_with_glia = annotations_mapped[inds1,]
var<-load(file=file.path(mappingFolder, paste(taxname, dataname, 'ann_map_roi_QC.Rdata', sep='_')))
annoNew_roi_with_glia = annoNew_roi_with_glia[,2:11]
colnames(annoNew_roi_with_glia) <- paste0(colnames(annoNew_roi_with_glia), '_wglia')
annotations_mapped2 <- merge(x = annoNew_sub, y = annoNew_roi_with_glia, by.x = "exp_component_name", by.y = 0, all=TRUE) 
annoNew_wglia <- annotations_mapped2[match(rownames(annoNew_roi_with_glia),annotations_mapped2$exp_component_name),]   # merge resorts things
save(annoNew_wglia, file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_roi_QC_wglia.Rdata', sep='_')))
write.csv(annoNew_wglia, file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_roi_QC_wglia.csv', sep='_')))


# Sampling counts
#load(file=file.path(mappingFolder,paste(taxname, dataname, 'sub_QC.Rdata"', sep='_')))
#compare to old full: var <- load('/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/NHP_BG_204_349_AIT115_ann_map_sub_QC.Rdata')


main_subclasses = c('STR D1 MSN GABA', 'STR D2 MSN GABA', 'STR Hybrid MSN GABA', 'STR Granular GABA', 'BG PVALB GABA', 'STR Cholinergic GABA', 'STR SST GABA', 'CN LHX6 LHX8 GABA', 'CN CGE LAMP5 GABA', 'CN MGE LAMP5 GABA')
# May now want to do this at the group level
#remaining_subclasses = setdiff(names(type_counts_Corr),main_subclasses) # TEMPORARILY SET TO CORR
#main_subclasses <- c(main_subclasses, remaining_subclasses)

dim(annoNew_sub)
y = annoNew_sub[[paste0(str_replace(subclass_colname,'_label',''),'_Corr')]]
type_counts_Corr = table(factor(y, levels=main_subclasses))
y = annoNew_sub[[paste0(str_replace(subclass_colname,'_label',''),'_Tree')]]
type_counts_Tree = table(factor(y, levels=main_subclasses))

#type_counts_Tree = table(annoNew_sub[paste0(str_replace(cluster_colname,'_label',''),'_Tree')])

# Plot NMS histogram for ROI
png(file.path(mappingFolder, paste(taxname, dataname, 'roi_NMS_hist.png', sep='_')), width = 500, height = 350)
hist(annoNew_roi_proj$marker_sum_norm_label)
dev.off()

#type_counts_Tree <- data.frame(type_counts_Tree) 
type_counts_Corr <- data.frame(type_counts_Corr) # TEMPORARILY SET TO CORR
#type_counts_Tree$Var1 <- factor(type_counts_Tree$Var1, levels = main_subclasses)
#type_counts_Tree$Var1 <- ordered(type_counts_Tree$Var1, levels = main_subclasses)
type_counts_Corr$Var1 <- ordered(type_counts_Corr$Var1, levels = main_subclasses) # TEMPORARILY SET TO CORR

colorkey <- read.csv(file = file.path('/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115', 'colortable.csv'))
#colors = list()
rownames(type_counts_Corr) = type_counts_Corr$Var1     # TEMPORARILY SET TO CORR
type_counts_Corr['color'] <- NA    # TEMPORARILY SET TO CORR
for (c in type_counts_Corr$Var1) {      # TEMPORARILY SET TO CORR
  tmp = colorkey['colz'][colorkey['subclass']==c]
  if (length(tmp)==0) {
    tmp = '#000000'
  } # Do we need this?
  #colors <- append (colors, tmp)
  type_counts_Corr[c,'color'] <- tmp       # TEMPORARILY SET TO CORR
  print(type_counts_Corr[c,'color'])        # TEMPORARILY SET TO CORR
  #print(colorkey$colz[colorkey['subclass']==c])
}
  
png(file.path(mappingFolder,'NHP_BG_AIT117_sampling_counts.png'), width = 1700, height = 1200)
tmp <- par("mar")
tmp[1] = tmp[1]+7
par(mar = tmp)
#barplot(names = type_counts_Tree$Var1, height = type_counts_Tree$Freq, las=2, col=type_counts_Tree$color)
ggplot(type_counts_Corr,aes(x = Var1, y = Freq)) +                     # TEMPORARILY SET TO CORR
       geom_bar(stat= 'Identity', fill = type_counts_Corr$color, alpha = 0.75) +
       #scale_fill_manual(values=type_counts_Tree$color) +
       xlab("") +
       ylab("Count") +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=40),
             theme(axis.text.x = element_text(angle = 55, vjust = 0.99, hjust=1, size=40),
             axis.text.y = element_text(size=40),
             axis.title=element_text(size=40, margin = margin(t = 0, r = 20, b = 0, l = 0)),
             #axis.line = element_line(size=1.5, color = 'midnightblue'),
             panel.background = element_blank()) +
       geom_hline(yintercept=10,linetype=2, size = 2, color = 'gray')
dev.off()

y = annoNew_sub$roi
y = gsub("_","",y)
y = gsub("([0-9])","",y)
y = gsub("STRdCP","STRd",y)
roi_counts = data.frame(table(y))
roi_counts = roi_counts[order(-roi_counts$Freq),]
roi_counts$y = factor(roi_counts$y, levels = roi_counts$y)

png(file.path(mappingFolder,'NHP_BG_AIT117_roi_counts.png'), width = 350, height = 400)
tmp <- par("mar")
tmp[1] = tmp[1]+7
par(mar = tmp)
ggplot(roi_counts, aes(x = y, y = Freq, fill = y)) +                     # TEMPORARILY SET TO CORR
  geom_bar(stat= 'Identity', alpha = 0.75) +
  xlab("") +
  ylab("Count") + guides(fill="none") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=24),
  theme(axis.text.x = element_text(angle = 55, vjust = 0.99, hjust=1, size=28),
        axis.text.y = element_text(size=24),
        axis.title=element_text(size=24, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        panel.background = element_blank())
dev.off()

# Striatal subclasses only (at least 5% of all cells are in dSTR or vSTR)
subclass = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
             "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
             "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
             "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
             "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
             "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")  

subclass = main_subclasses

# cumsum
df_NMS <- subset(annoNew_roi_proj, Subclass_Corr %in% c('SN_STH', 'SLC17A7-SATB2', 'MEIS2', 'D1-Matrix', 'D2-Matrix', 'PVALB-COL19A1-ST18'))
png(file.path(mappingFolder,'NHP_BG_AIT117_cumNMS_MEIS2_SubclassCorr_scoreCorr_roi_proj.png'), width = 500, height = 250)
#ggplot(df_NMS, aes(marker_sum_norm_label, color=Subclass_Tree)) + 
ggplot(df_NMS, aes(score.Corr, color=Subclass_Corr)) + 
  stat_ecdf(geom = "step", size=1) + 
  labs(x="score.Corr", y="Cumulative Fraction") +
  scale_y_continuous(breaks=seq(0,1,0.1), labels = seq(0,100,10)) + 
  #scale_color_manual(values=c("#96ceb4", "#ff6f69", "#ffcc5c", "#90697c"), name='Subclass') + 
  scale_color_manual(values=c("#758f0b", "#ff6f69", "#ffcc5c", "#90697c", "#007c8f", "#cf0690"), name='Subclass') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
dev.off()

df_NMS <- subset(annoNew_sub, Subclass_Corr %in% c('MEIS2', 'LHX6-TAC3-PLPP4', 'TAC3-LHX8-PLPP4', 'D1-ICj', 'D1-NUDAP', 'D1-ShellOT', 'D2-ShellOT'))
df_NMS$Region2 = gsub("_", "", df_NMS$Region)
png(file.path(mappingFolder,'NHP_BG_AIT117_cumNMS_MEIS2_SubclassCorr_sub_roi_dist.png'), width = 500, height = 350)
#ggplot(df_NMS, aes(marker_sum_norm_label, color=Subclass_Tree)) + 
ggplot(df_NMS, aes(fill=Region2, x=Subclass_Corr)) + 
  geom_bar(position="stack", stat="count") +
  labs(x="", y="count") +
  #scale_color_manual(values=c("#ff6f69", "#ffcc5c", "#90697c", "#007c8f", "#cf0690"), name='Subclass') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12, angle=90),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=16, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
dev.off()

df_NMS <- subset(annoNew_sub, Subclass_Tree %in% c('MEIS2', 'LHX6-TAC3-PLPP4', 'TAC3-LHX8-PLPP4', 'D1-ICj', 'D1-NUDAP', 'D1-ShellOT', 'D2-ShellOT'))
df_NMS$Region2 = gsub("_", "", df_NMS$Region)
png(file.path(mappingFolder,'NHP_BG_AIT117_cumNMS_MEIS2_SubclassTree_sub_roi_dist.png'), width = 500, height = 350)
#ggplot(df_NMS, aes(marker_sum_norm_label, color=Subclass_Tree)) + 
ggplot(df_NMS, aes(fill=Region2, x=Subclass_Tree)) + 
  geom_bar(position="stack", stat="count") +
  labs(x="", y="count") +
  #scale_color_manual(values=c("#ff6f69", "#ffcc5c", "#90697c", "#007c8f", "#cf0690"), name='Subclass') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12, angle=90),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=16, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
dev.off()

postpatch1 = annoNew_sub['postPatch_classification'][annoNew_sub['level3.subclass_Tree'] == 'SN_STH']
table(postpatch1)/length(postpatch1)
postpatch2 = annoNew_sub['postPatch_classification'][annoNew_sub['level3.subclass_Tree'] == 'D1-Matrix']
table(postpatch2)/length(postpatch2)

# Count cells with ephys
df_ephys = read.csv(file=file.path("NHP_ephys_features_20240822.csv"))
df_id = read.csv("custom_report_20240429.csv")

df2 = merge(annoNew_sub, df_id, by.x='cell_name', by.y='cell_specimen_name.', all.x = FALSE, all.y = FALSE)
# checked these numbers against python by merging against annoNew_roi instead
df3 = merge(df2, df_ephys, by.x ='cell_specimen_id.', right_on='cell_name', all.x = FALSE, all.y = FALSE)
sum(is.na(df3['tau']))
sum(is.na(df3['adapt_hero']))
sum(is.na(df3['upstroke_downstroke_ratio_short_square']))
grep("sag", colnames(df3))
df_ephys = df3[265:358]
df_ephys <- df_ephys[,colSums(is.na(df_ephys))<nrow(df_ephys)]  # drop any rows with all Na; nothing dropped

df_short = subset(type_counts_Corr, select = -c(color))   # TEMPORARILY SET TO CORR
df_short['source'] = 'Tx'
#type_counts_ephys = table(df3$level3.subclass_Corr)   # TEMPORARILY SET TO CORR
type_counts_ephys = table(factor(df3$Subclass_Corr, levels=main_subclasses)) 
type_counts_ephys <- data.frame(type_counts_ephys)
type_counts_ephys['source'] = 'Ephys'
rownames(type_counts_ephys) = type_counts_ephys$Var1
type_counts = rbind(df_short,type_counts_ephys)
type_counts$source <- ordered(type_counts$source, levels = c('Tx','Ephys'))

png(file.path(mappingFolder,'NHP_BG_AIT117_sampling_counts_wephys.png'), width = 1800, height = 1200)
tmp <- par("mar")
tmp[1] = tmp[1]+7
par(mar = tmp)
#barplot(names = type_counts_Tree$Var1, height = type_counts_Tree$Freq, las=2, col=type_counts_Tree$color)
ggplot(type_counts,aes(x= Var1, y = Freq, fill = source)) +
  geom_bar(stat= 'Identity',position="dodge") +
  #scale_fill_manual(values=type_counts_Tree$color) +
  xlab("") +
  ylab("Count") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 36), 
        axis.text.y = element_text(size=36),
        axis.title=element_text(size=40, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text=element_text(size=40),
        legend.title=element_text(size=40)) +
  geom_hline(yintercept=10,linetype=2, size = 2, color = 'gray') +
  scale_fill_manual(values = c("#62b7c4", "#d455ac"))
dev.off()

unique(annotations_mapped$roi)

anno_mapped_roi = annoNew[inds1,]
query.data_roi = query.data[,inds1]
write.csv(anno_mapped_roi, file.path(mappingFolder, paste(taxname, dataname, 'ann_map_roi.csv', sep='_')), row.names=FALSE)
type_counts_Corr_roi = table(anno_mapped_roi$level3.subclass_Corr)
type_counts_Tree_roi = table(anno_mapped_roi$level3.subclass_Tree)

#anno_mapped_sub = annotations_mapped[inds1&inds3&inds4&inds6,]
query.data_sub = query.data[,inds1&inds3&inds4&inds6]

rownames(annoNew) <- annoNew$exp_component_name   # After running apply_PatchseqQC
anno_mapped_sub <- annoNew[inds1&inds3&inds4&inds6,]  

type_counts_Corr_QC = table(anno_mapped_sub$level3.subclass_Corr)
type_counts_Tree_QC = table(anno_mapped_sub$level3.subclass_Tree)
write.csv(anno_mapped_sub, file.path(mappingFolder, paste(taxname, dataname, 'ann_map_sub_QC.csv', sep='_')), row.names=FALSE)

# QC files for Rachel
# Optional load annoNew from another run
# Run on server:
load(file=file.path(mappingFolder,"NHP_BG_AIT_117_RSC-204-373_roi_QC.Rdata"))
annoNew <- annoNew_roi
# Run on laptop:
#load(file="/Users/xiaoping.liu/celltypes/NHP_BG_anal/NHP_BG_AIT_115/204_359/NHP_BG_204_359_AIT115_ann_map_full_QC.Rdata")
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

desired_columns = c('exp_component_name', 'cell_name', 'cell_id', 'Subclass_Corr', 
                    'Subclass_Tree', 'rna_amplification_pass_fail', 'compound_qc_pass', 
                    'BG_ROI', 'roi', 'species', 'postPatch_classification', 'acute', 'Virus', 
                    'creCell', 'revisit')
# Or striatal ROI?
anno_morpho = annoNew[desired_columns]
write.csv(anno_morpho, file.path(mappingFolder,"NHP_BG_204_373_AIT117_anno_morpho.csv"))


Tax_Ca_Pu = AIT.anndata$obs[AIT.anndata$obs$roi_label %in% c('Macaque CaB', 'Macaque CaH', 'Macaque CaT', 'Macaque PuC', 
                'Macaque PuPV', 'Macaque PuR'),c('level1.class_label', 'level3.subclass_label')]
table(Tax_Ca_Pu$level1.class_label)
table(Tax_Ca_Pu$level3.subclass_label)

Ps_Ca_Pu = annoNew_sub[annoNew_sub$roi %in% c('STRdPu', 'STRdCa', 'STRd', 'STRdCP', 'STRd_CP', 'STRd_Pu', 'STRd_Ca'),]
table(Ps_Ca_Pu$level1.class_Corr)
table(Ps_Ca_Pu$level3.subclass_Corr)


# Ephys feature analysis in R
efeats <- read.csv('/home/xiaoping.liu/Desktop/NHP_ephys_features_20240430.csv')
ids <- read.csv('/home/xiaoping.liu/Desktop/custom_report_20240610.csv')
df <- merge(efeats, ids, by.x="cell_name", by.y="cell_specimen_id.")

load(paste0(data_dir, paste0(data_fn, "_cpm.Rdata")))
load(paste0(data_dir, paste0(data_fn, "_samp.dat.Rdata")))

query.metadata <- samp.dat
counts      <- cpmR   # Genes are rows, samples are columns

query.counts   <- counts
query.data   <- logCPM(query.counts)

df <- merge(df, query.metadata, by.x= 'cell_specimen_name.', by.y='cell_name')
df <- merge(df, query.data['PVALB',], by.x="exp_component_name", by.y='row.names')

# Subset to QC sub
df <- merge(df, annoNew_sub, by = 'exp_component_name')

# Subset to PTHLH neurons
df <- df[df$Subclass_Corr == "PVALB-COL19A1-ST18",]

# Calculate correlation
result = cor.test(df$y, df$width_rheo, method = "kendall", use="pairwise.complete.obs")
result = cor.test(df$y, df$mean_isi_hero, method = "kendall", use="pairwise.complete.obs")

png(file.path(mappingFolder, 'width_rheo_vs_PVALB_sub_nolog.png'), width = 500, height = 350)
ggplot(df, aes(x=y, y=width_rheo)) + 
  geom_point() + xlab('CPM_PVALB')
dev.off()

png(file.path(mappingFolder, 'mean_isi_hero_vs_PVALB_nolog.png'), width = 500, height = 350)
ggplot(df, aes(x=y, y=mean_isi_hero)) + 
#  geom_point() + xlab('log2(CPM_PVALB+1)')
  geom_point() + xlab('CPM_PVALB')
dev.off()


