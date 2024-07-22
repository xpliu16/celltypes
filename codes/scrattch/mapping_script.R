suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)
library (dplyr)
library(stringr)

#START

source("/home/xiaoping.liu/scrattch/mapping/run_mappings.R") # on cluster

# CrossAreal taxonomies

#region = "MTG"
#roi_strs = "TCx|Tcx|TEa"

#region = "V1"
#roi_strs = "OCx|VISp"

#region = "M1"
#roi_strs = "FCx|MOp"

region = "DLPFC"
roi_strs = "FCx|dlPFC"

data_fn = "20240321_RSC-122-359_human_patchseq_star2.7"
  
run_mappings(refFolder = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_",region),
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/Human_Cortical", 
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/",
             data_fn = data_fn,
             mode = paste0("CrossAreal_",region,"_patchseq"),
             h5ad_fn = paste0("CrossAreal_",region,".h5ad"), 
             class_colname = "Class_label",
             neigh_colname = NULL,
             subclass_colname = "CrossArea_subclass_label", 
             cluster_colname = "CrossArea_cluster_label", 
             proj_strs = "hIVSCC-MET",
             roi_strs = roi_strs,
             off_target = "non-neuronal"
)
# For reloading data and exploratory work:
refFolder = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_",region)
mappingFolder = "/home/xiaoping.liu/scrattch/mapping/Human_Cortical"
data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/"
mode = paste0("CrossAreal_",region,"_patchseq")
h5ad_fn = paste0("CrossAreal_",region,".h5ad")
class_colname = "Class_label"
neigh_colname = NULL
subclass_colname = "CrossArea_subclass_label" 
cluster_colname = "CrossArea_cluster_label"
proj_strs = "hIVSCC-MET"
off_target = "non-neuronal"

# BG
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115_NCBI"
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115_NCBI"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/"
mode = 'patchseq'
data_fn = "20240425_RSC-204-361_macaque_patchseq_star2.7"
h5ad_fn = NULL  # AI_taxonomy.h5ad 
class_colname = 'level1.class_label'
neigh_colname = 'level2.neighborhood_label'
subclass_colname = 'level3.subclass_label'
cluster_colname = 'cluster_label'
proj_strs = "qIVSCC-MET"
roi_strs = "STR"
off_target = 'NN'

run_mappings(refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_116",
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116", 
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/",
             data_fn = "20240621_RSC-204-366_macaque_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = NULL, 
             class_colname = 'Class_label',
             neigh_colname = 'Neighborhood_label',
             subclass_colname = 'Subclass_label', 
             cluster_colname = 'Cluster_label', 
             proj_strs = "qIVSCC-MET",
             roi_strs = "STR",
             off_target = "NN"
)
refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_116" 
mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116"  
data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/" 
data_fn = "20240520_RSC-204-363_macaque_patchseq_star2.7" 
mode = 'patchseq'                                                                                  
h5ad_fn = NULL  
class_colname = 'Class_label' 
neigh_colname = 'Neighborhood_label' 
subclass_colname = 'Subclass_label'  
cluster_colname = 'Cluster_label'  
proj_strs = "qIVSCC-MET" 
roi_strs = "STR"
off_target = "NN"


run_mappings(refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115_NCBI",
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115_NCBI", 
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/",
             data_fn = "20240425_RSC-204-361_macaque_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = NULL, 
             class_colname = 'level1.class_label',
             neigh_colname = 'level2.neighborhood_label',
             subclass_colname = 'level3.subclass_label', 
             cluster_colname = 'cluster_label', 
             proj_strs = "qIVSCC-MET",
             roi_strs = "STR",
             off_target = "NN"
)

# MTG taxonomy

run_mappings(refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Macaque_NCBI",
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/GreatApes_Macaque_NCBI", 
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/",
             data_fn = "20240425_RSC-204-361_macaque_patchseq_star2.7",
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

run_mappings(refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Human",
             mappingFolder = "/home/xiaoping.liu/scrattch/mapping/GreatApes_Human", 
             data_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/",
             data_fn = "20240606_RSC-122-365_human_patchseq_star2.7",
             mode = 'patchseq',
             h5ad_fn = 'GreatApes_Macaque_NCBI.h5ad', 
             class_colname = 'class_label',
             neigh_colname = 'neighborhood_label',
             subclass_colname = 'subclass_label', 
             cluster_colname = 'cluster_label', 
             proj_strs = "qIVSCC-MET",
             roi_strs = "",
             off_target = "glia"
)



run_mappings <- function(refFolder, mappingFolder, data_dir, data_fn, mode, 
                         h5ad_fn = NULL, class_colname, neigh_colname, 
                         subclass_colname, cluster_colname, proj_strs, roi_strs,
                         off_target) {
    # This remakes taxonomy from feather files and reloads, very slow, use read_h5ad instead for the time being
    if (is.null(h5ad_fn)) {
        AIT.anndata <- loadTaxonomy(refFolder)
    } else {
        AIT.anndata <- loadTaxonomy(refFolder, h5ad_fn)
    }

    tryCatch({
      AIT.anndata = mappingMode(AIT.anndata, mode=mode)
      }, error = function(err) {
        ## Add in the off.target annotation.
        AIT.anndata$obs['off_target'] = AIT.anndata$obs[class_colname]
        
        # Setup the taxonomy for patchseqQC to infer off.target contamination
        AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                            mode.name = mode, ## Give a name to off.target filterd taxonomy
                                            subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation
                                            subclass.column = subclass_colname, ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                            class.column = class_colname, ## The column by which off-target types are determined
                                            off.target.types = off_target, ## The off-target class.column labels for patchseqQC
                                            num.markers = 50, ## Number of markers for each annotation in `class_label`
                                            taxonomyDir = refFolder)
        AIT.anndata = mappingMode(AIT.anndata, mode=mode)
      })

    load(paste0(data_dir, paste0(data_fn, "_cpm.Rdata")))
    load(paste0(data_dir, paste0(data_fn, "_samp.dat.Rdata")))

    query.metadata <- samp.dat
    counts      <- cpmR   # Genes are rows, samples are columns

    query.counts   <- counts
    query.data   <- logCPM(query.counts)

    # Put annotations and counts in the same order
    query.metadata <- query.metadata[match(colnames(query.data),query.metadata$exp_component_name),] 
    rownames(query.metadata) <- query.metadata$exp_component_name  

    # If necessary, trim data down to only unique cellnames
    
    ## Get marker genes from dendrogram
    #addDendrogramMarkers()  # When running for first time
    #load(AIT.anndata$uns$dend)
    #dend = reference$dend
    #dend <- readRDS(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
    dend <- json_to_dend(fromJSON(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]]))
    #dend <- readRDS(file.path(refFolder, mode, "dend.RData"))  # Or whatever mode is

    #dend_file = paste0(refFolder,'/dend.RData')
    #dend <- readRDS(dend_file)
    allMarkers = unique(unlist(get_dend_markers(dend)))

    ## Alternative get genes from MERFISH panel
    #panel_df = read.csv(file.path(panelFolder, "AIT_115_MERFISH_panel.xlsx"))

    ## Subset query data to just those markers
    query.data = query.data[intersect(rownames(query.data), allMarkers),]

    ## Identify the offtarget cell types manually.
    print(unique(AIT.anndata$obs[class_colname]))

    query.mapping_obj <- taxonomy_mapping(AIT.anndata= AIT.anndata,
                                    query.data = query.data, 
                                    corr.map   = TRUE, # Flags for which mapping algorithms to run
                                    tree.map   = TRUE, 
                                    seurat.map = FALSE, 
                                    label.cols = c(class_colname, neigh_colname, subclass_colname, cluster_colname)
    )
    a <- strsplit(refFolder,'/')[[1]]
    taxname <- a[length(a)]
    b <- strsplit(data_fn, '_')[[1]]
    dataname <- b[2]
    
    if (! file.exists(mappingFolder)) {
      dir.create(mappingFolder)
    }

    save(query.mapping_obj, file=file.path(mappingFolder, paste(taxname, dataname, 'mapping.Rdata', sep='_')))
    query.mapping = getMappingResults(query.mapping_obj)
    write.csv(query.mapping, file.path(mappingFolder, paste(taxname, dataname, 'mapping.csv', sep='_')), row.names=FALSE)
    
    # Variable renaming
    #clusters  <- unique(query.mapping$cluster)   
    clusters <- unique(AIT.anndata$uns$clusterInfo[cluster_colname])
    #subclass_levels <- unique(query.mapping$subclass_Corr)
    #subclass_levels <- unique(AIT.anndata$uns$clusterInfo$subclass_label)

    annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = 0, all=TRUE) 
    #annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = "exp_component_name", all=TRUE) 
    # Alt hack:
    #rownames(query.mapping) = rownames(query.metadata)
    #annotations_mapped <- merge(x = query.mapping, y = query.metad ata, by.x = 0, by.y = 0, all=TRUE) 
    # annotations_mapped$clusters <- factor(annotations_mapped$cluster_Tree, levels=clusters)
    # annotations_mapped$class   <- annotations_mapped$class_Tree
    # annotations_mapped$class[!is.element(annotations_mapped$class,c("IN","MSN"))] = "Non-neuronal"
    annotations_mapped <- annotations_mapped[match(rownames(query.metadata),annotations_mapped$exp_component_name),]   # merge resorts things
    rownames(annotations_mapped) <- annotations_mapped$exp_component_name
    #type_counts_Corr = table(annotations_mapped$supertype_Corr)
    #type_counts_Tree = table(annotations_mapped$supertype_Tree)
    type_counts_Corr = table(annotations_mapped$level3.subclass_Corr)
    type_counts_Tree = table(annotations_mapped$level3.subclass_Tree)

    write.csv(annotations_mapped, file=file.path(mappingFolder, paste(taxname, dataname, 'ann_map_full.csv', sep='_')), row.names=FALSE)
    save(annotations_mapped, type_counts_Corr, type_counts_Tree, file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_full.Rdata', sep='_')))

    buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                        mappingFolder  = file.path(mappingFolder, paste(taxname, dataname, 'map_full', sep='_')),
                        query.data     = counts,  # Don't need log-normalized data here
                        query.metadata = query.metadata,
                        query.mapping  = query.mapping_obj,
                        doPatchseqQC   = TRUE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
    )

    dir(file.path(mappingFolder, paste(taxname, dataname, 'map_full', sep='_')))


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
                            query.metadata = annotations_mapped, ## Results of the previous mapping or AIT.anndata$obs, no mapping is required.
                            verbose=FALSE)
    # Ran the contents of the function directly in R to avoid error in /R/patchseq_output.R
    save(annoNew, file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_full_QC.Rdata', sep='_')))
    write.csv(annoNew, file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_full_QC.csv', sep='_')))
    rownames(annoNew) = annoNew$exp_component_name
    
    # annotations_mapped$cluster <- factor(annotations_mapped$clusters, levels=clusters)  # Make into discrete levels
    # see all the regions: unique(sapply(annoNew$roi, function(x) strsplit(x, '_')[[1]][1]))
    #annoNew_hIVSCC = annoNew[grepl("hIVSCC-MET", annoNew$cell_specimen_project),]
    #inds1 = grepl("STR",annoNew$roi)
    inds0 = grepl(proj_strs, annoNew$cell_specimen_project)
    inds1 = grepl(roi_strs, annoNew$roi)
    #inds1 = ifelse(grepl("STR|PALGPi|HYSTN",annoNew$roi), TRUE,FALSE)
    #inds2 = annotations_mapped$library_prep_pass_fail == "Pass"   # Chucks good samples
    inds3 = annoNew$Genes.Detected >= 1000
    inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
    #inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
    #inds6 = annotations_mapped$score.Corr > 0.6
    inds6 = annoNew$marker_sum_norm_label >= 0.6
    
    annoNew_roi = annoNew[inds1,]
    save(annoNew_roi, file=file.path(mappingFolder, paste(taxname, dataname, 'roi_QC.Rdata', sep='_')))
    write.csv(annoNew_roi, file=file.path(mappingFolder, paste(taxname, dataname, 'roi_QC.csv', sep='_')))
    
    annoNew_roi_proj = annoNew[inds0&inds1,]
    save(annoNew_roi_proj, file=file.path(mappingFolder, paste(taxname, dataname, 'roi_proj_QC.Rdata', sep='_')))
    
    annoNew_sub = annoNew[inds0&inds1&inds3&inds4&inds6,]
    save(annoNew_sub, file=file.path(mappingFolder,paste(taxname, dataname, 'sub_QC.Rdata', sep='_')))
    write.csv(annoNew_sub, file=file.path(mappingFolder,paste(taxname, dataname, 'sub_QC.csv', sep='_')))
}

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


main_subclasses = c('D1-Matrix', 'D2-Matrix', 'D1-Striosome', 'D2-Striosome', 'D1-ShellOT', 'D2-ShellOT', 'D1D2-Hybrid', 'D2-Hybrid-MCHR2', 'D1-NUDAP', 'D1-ICj', 'PVALB-COL19A1-ST18', 'SST_Chodl', 'CCK-FBXL7', 'CHAT', 'CCK-VIP-TAC3', 'LHX6-TAC3-PLPP4', 'TAC3-LHX8-PLPP4', 'MEIS2')
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

colorkey <- read.csv(file = file.path(mappingFolder, 'colortable.csv'))
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

png(file.path(mappingFolder,'NHP_BG_AIT115_sampling_counts.png'), width = 2000, height = 1200)
tmp <- par("mar")
tmp[1] = tmp[1]+7
par(mar = tmp)
#barplot(names = type_counts_Tree$Var1, height = type_counts_Tree$Freq, las=2, col=type_counts_Tree$color)
ggplot(type_counts_Corr,aes(x = Var1, y = Freq)) +                     # TEMPORARILY SET TO CORR
       geom_bar(stat= 'Identity', fill = type_counts_Corr$color, alpha = 0.75) +
       #scale_fill_manual(values=type_counts_Tree$color) +
       xlab("") +
       ylab("Count") +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=40),
             axis.text.y = element_text(size=40),
             axis.title=element_text(size=40, margin = margin(t = 0, r = 20, b = 0, l = 0)),
             #axis.line = element_line(size=1.5, color = 'midnightblue'),
             panel.background = element_blank()) +
       geom_hline(yintercept=10,linetype=2, size = 2, color = 'gray')
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
png(file.path(mappingFolder,'NHP_BG_AIT116_cumNMS_MEIS2_SubclassCorr_scoreCorr_roi_proj.png'), width = 500, height = 250)
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
png(file.path(mappingFolder,'NHP_BG_AIT116_cumNMS_MEIS2_SubclassCorr_sub_roi_dist.png'), width = 500, height = 350)
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
png(file.path(mappingFolder,'NHP_BG_AIT116_cumNMS_MEIS2_SubclassTree_sub_roi_dist.png'), width = 500, height = 350)
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
df_ephys = read.csv(file=file.path("NHP_ephys_features_20240430.csv"))
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
type_counts_ephys = table(factor(df3$level3.subclass_Corr, levels=main_subclasses)) 
type_counts_ephys <- data.frame(type_counts_ephys)
type_counts_ephys['source'] = 'Ephys'
rownames(type_counts_ephys) = type_counts_ephys$Var1
type_counts = rbind(df_short,type_counts_ephys)
type_counts$source <- ordered(type_counts$source, levels = c('Tx','Ephys'))

png(file.path(mappingFolder,'NHP_BG_AIT115_sampling_counts_wephys.png'), width = 1800, height = 1200)
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


dir.create(file.path(mappingFolder, paste(taxname, dataname, 'map_sub_patchseqQC', sep='_')))
#unlink(file.path(mappingFolder, "NHP_BG_204_358_AIT115_map_sub_patchseqQC/dend.RData"))

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = file.path(mappingFolder, paste(taxname, dataname, 'map_sub_patchseqQC', sep='_')),
                      query.data     = query.counts_sub,  # Don't need log-normalized data here
                      query.metadata = query.metadata_sub,
                      query.mapping  = query.mapping_sub,
                      doPatchseqQC   = TRUE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

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
load(file=file.path(mappingFolder,"NHP_BG_AIT_115_NCBI_RSC-204-359_roi_QC.Rdata"))
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
annoNew$revisit1 = (annoNew$rna_amplification_pass_fail=="Fail") & (annoNew$compound_qc_pass == TRUE) 

desired_columns = c('exp_component_name', 'cell_name', 'cell_id', 'level3.subclass_Corr', 
                    'level3.subclass_Tree', 'rna_amplification_pass_fail', 'compound_qc_pass', 
                    'BG_ROI', 'roi', 'species', 'postPatch_classification', 'acute', 'Virus', 
                    'creCell')
# Or striatal ROI?
anno_morpho = annoNew[desired_columns]
write.csv(anno_morpho, file.path(mappingFolder,"NHP_BG_204_359_AIT115_anno_morpho.csv"))


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

png(file.path(mappingFolder, 'width_rheo_vs_PVALB.png'), width = 500, height = 350)
ggplot(df, aes(x=y, y=width_rheo)) + 
  geom_point() + xlab('log2(CPM_PVALB+1)')
dev.off()

png(file.path(mappingFolder, 'mean_isi_hero_vs_PVALB.png'), width = 500, height = 350)
ggplot(df, aes(x=y, y=mean_isi_hero)) + 
  geom_point() + xlab('log2(CPM_PVALB+1)')
dev.off()


