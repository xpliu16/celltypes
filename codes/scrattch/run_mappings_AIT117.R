suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.taxonomy)
library(scrattch.mapping)
library(scrattch.patchseq)
library (dplyr)
library(stringr)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")
#reticulate::use_python("/usr/bin/python3")

run_mappings <- function(refFolder, mappingFolder, data_dir, data_fn, mode, 
                         h5ad_fn = NULL, hierarchy, proj_strs, roi_strs, 
                         off_target, off_target_level) {
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
    #AIT.anndata$obs['off_target'] = AIT.anndata$obs[off_target_level]
    
    # Setup the taxonomy for patchseqQC to infer off.target contamination
    AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                        mode.name = mode, ## Give a name to off.target filterd taxonomy
                                        subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation
                                        subclass.column = subclass_colname, ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                        class.column = class_colname, ## The column by which off-target types are determined
                                        off.target.types = off_target, ## The off-target class.column labels for patchseqQC
                                        subclass.subsample = 100,
                                        num.markers = 50, ## Number of markers for each annotation in `class_label`
                                        taxonomyDir = refFolder)
    AIT.anndata = mappingMode(AIT.anndata, mode=mode)
  })
  
  load(paste0(data_dir, paste0(data_fn, "_cpm.Rdata")))
  load(paste0(data_dir, paste0(data_fn, "_samp.dat.Rdata")))
  
  #Check for uniqueness
  if (length(unique(colnames(cpmR)))!=dim(cpmR)[2]){
    print("Error: sample names are not unique")
    return()
  } 
  
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
  dend <- json_to_dend(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
  #dend <- readRDS(file.path(refFolder, mode, "dend.RData"))  # Or whatever mode is
  
  #dend_file = paste0(refFolder,'/dend.RData')
  #dend <- readRDS(dend_file)
  allMarkers = unique(unlist(get_dend_markers(dend)))

  # Filter out "trouble" genes (but beware tree mapping may not work)
  excl1 = read.csv("/home/xiaoping.liu/Jeremy_exclude_gene_lists/sex_genes.txt")
  #excl2 = read.csv("/home/xiaoping.liu/Jeremy_exclude_gene_lists/mito_genes.txt")
  #excl3 = read.csv("/home/xiaoping.liu/Jeremy_exclude_gene_lists/activity_genes.txt")
  #excl4 = read.csv("/home/xiaoping.liu/Jeremy_exclude_gene_lists/is_immune_pvalb.txt")
  #foo<-unique(unlist(unlist(list(excl1, excl2, excl3, excl4))))
  excl<-unique(unlist(unlist(list(excl1))))
  allMarkers = setdiff(allMarkers, excl) # 1382 to 1347 (sex genes) 1313 (all trouble genes)

  ## Alternative get genes from MERFISH panel
  #panel_df = read.csv(file.path(panelFolder, "AIT_115_MERFISH_panel.xlsx"))
  
  ## Subset query data to just those markers
  query.data = query.data[intersect(rownames(query.data), allMarkers),]

   ## Compute top 1000 binary marker genes for clusters (or use a pre-existing vector)
  keep.cells   = !AIT.anndata$uns$filter[[mode]]
  unique(AIT.anndata$obs$Subclass_label[keep.cells])    # check off-target filtering
  binary.genes = top_binary_genes(t(AIT.anndata$X[AIT.anndata$obs_names[keep.cells],]), AIT.anndata$obs$cluster_label[keep.cells], 1500)
  # OR as.matrix(AIT.anndata$X)
  # cache binary genes to speed up and examine
  # write.csv(binary.genes, "/home/xiaoping.liu/scrattch/mapping/binarygenes.csv", row.names=FALSE)

  ## Update the anndata object with this gene set
  AIT.anndata  = updateHighlyVariableGenes(AIT.anndata,binary.genes)
  #AIT.anndata = updateHighlyVariableGenes(AIT.anndata, allMarkers)
  # Check
  sum(AIT.anndata$var$highly_variable_genes_patchseq)

  
  ## Identify the offtarget cell types manually.
  print(unique(AIT.anndata$obs[class_colname]))
  
  query.mapping_obj <- taxonomy_mapping(AIT.anndata= AIT.anndata,
                                        query.data = query.data, 
                                        label.cols=hierarchy,
                                        corr.map=TRUE,
                                        hierarchical.map=FALSE,
                                        tree.map=TRUE,
                                        seurat.map=TRUE
  )
  
  a <- strsplit(refFolder,'/')[[1]]
  taxname <- a[length(a)]
  b <- strsplit(data_fn, '_')[[1]]
  dataname <- b[2]

  # Dataname modifier for testing mapping strategy variants
  mappingFolder <- file.path(mappingFolder,"newdocker_binary1500_noexclusions_20250219")
  dir.create(mappingFolder, showWarnings = FALSE)
  
  save(query.mapping_obj, file=file.path(mappingFolder, paste(taxname, dataname, 'mapping.Rdata', sep='_')))
  query.mapping = getMappingResults(query.mapping_obj, scores = TRUE)
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
  #type_counts_Corr = table(annotations_mapped$level3.subclass_Corr)
  #type_counts_Tree = table(annotations_mapped$level3.subclass_Tree)
  type_counts_Corr = table(annotations_mapped$Group_Corr)
  type_counts_Tree = table(annotations_mapped$Group_Tree)
  type_counts_Seurat = table(annotations_mapped$Group_Seurat)

  write.csv(annotations_mapped, file=file.path(mappingFolder, paste(taxname, dataname, 'ann_map_full.csv', sep='_')), row.names=FALSE)
  save(annotations_mapped, type_counts_Corr, type_counts_Tree, file=file.path(mappingFolder,paste(taxname, dataname, 'ann_map_full.Rdata', sep='_')))

  createShiny(AIT.anndata,
              shinyDir = file.path(mappingFolder, paste(taxname, dataname, 'map_full', sep='_')), # Replace location with desired location for shiny directory output
              metadata_names = NULL)
  
  #buildMappingDirectory(AIT.anndata    = AIT.anndata, 
  #                    mappingFolder  = file.path(mappingFolder, paste(taxname, dataname, 'map_full', sep='_')),
  #                    query.data     = counts,  # Don't need log-normalized data here
  #                    query.metadata = query.metadata,
  #                    query.mapping  = query.mapping_obj,
  #                    doPatchseqQC   = TRUE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
  #                    return.metrics = TRUE
  #)

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
  #inds3 = annoNew$Genes.Detected >= 1000
  inds3 = annoNew$Genes.Detected >= 1500
  inds4 = annoNew$percent_reads_aligned_total >= 25      # Very conservative, but looks like nothing chucked improperly on UMAP
  #inds5 = annotations_mapped$percent_reads_aligned_to_introns > 25   # Chucks good samples
  #inds6 = annotations_mapped$score.Corr > 0.6
  #inds6 = annoNew$marker_sum_norm_label >= 0.55
  #inds6 = annoNew$score.Tree >= 0.4
  inds6 = annoNew$marker_sum_norm_label >= 0.6
  #inds6 = annoNew$score.Corr>=0.4
  
  cols <- colnames(annoNew)
  first_cols <- c("Row.names", "cell_name", "score.Corr", "score.Tree", "score.Seurat", 
                  "Group_Corr", "Cluster_Corr",
                  "Group_Tree",	"Cluster_Tree",
                  "Group_Seurat", "Cluster_Seurat",
                  "exp_component_name",	"marker_sum_norm_label", "Genes.Detected",
                  "percent_reads_aligned_to_introns",	"percent_reads_aligned_total",	
                  "postPatch", "endPipetteR", "Virus", "creCell",	"Region",	
                  "postPatch_classification",	"go_no_go_63x",	"image_series_63x_qc",
                  "rna_amplification_pass_fail", "percent_cdna_longer_than_400bp",
                  "amplified_quantity_ng")
  reordered_cols = c(first_cols, setdiff(cols, first_cols))
  annoNew = annoNew[,reordered_cols]
  annoNew_roi = annoNew[inds1,]
  save(annoNew_roi, file=file.path(mappingFolder, paste(taxname, dataname, 'roi_QC.Rdata', sep='_')))
  write.csv(annoNew_roi, file=file.path(mappingFolder, paste(taxname, dataname, 'roi_QC.csv', sep='_')))
  
  annoNew_roi_proj = annoNew[inds0&inds1,]
  save(annoNew_roi_proj, file=file.path(mappingFolder, paste(taxname, dataname, 'roi_proj_QC.Rdata', sep='_')))
  
  annoNew_sub = annoNew[inds0&inds1&inds3&inds4&inds6,]
  save(annoNew_sub, file=file.path(mappingFolder,paste(taxname, dataname, 'sub_QC_NMS06.Rdata', sep='_')))
  write.csv(annoNew_sub, file=file.path(mappingFolder,paste(taxname, dataname, 'sub_QC_NMS06.csv', sep='_')))
  
  type_counts_Corr_sub = table(annoNew_sub$Group_Corr)
  type_counts_Tree_sub = table(annoNew_sub$Group_Tree)
  type_counts_Seurat_sub = table(annoNew_sub$Group_Seurat)
}

# If reloading
#rownames(annoNew_sub) <- annoNew_sub[,'exp_component_name']

query.data[c('SST','NPY', 'NOS1', 'TACR1', 'SCN9A'),rownames(annoNew_sub)[annoNew_sub$Group_Corr=='STR SST Chodl GABA']]
query.data['SST',rownames(annoNew_sub)[annoNew_sub$Group_Tree=='STR SST Chodl GABA']]
query.data['SST',rownames(annoNew_sub)[annoNew_sub$Group_Tree=='STR SST RSPO2 GABA']]
query.data[c('CPNE4', 'KCNT2', 'DRD3', 'NTN1', 'PPP1R1B') ,rownames(annoNew_sub)[annoNew_sub$Group_Tree=='OT D1-ICj']]
query.data[c('CHST9', 'RXFP1', 'B3GAT2', 'PPP1R1B', 'DRD2') ,rownames(annoNew_sub)[annoNew_sub$Group_Tree=='STRv D1 NUDAP']]
query.data[c('CHST9', 'RXFP1', 'B3GAT2', 'PPP1R1B', 'DRD2') ,rownames(annoNew_sub)[annoNew_sub$Group_Tree=='STRd D1D2 Hybrid']]
query.data[c('DRD1','KCNIP1', 'STXBP6', 'DRD2') ,rownames(annoNew_sub)[annoNew_sub$Group_Tree=='STRd D1 Striosome']]
query.data[c('DRD1','KCNIP1', 'STXBP6', 'DRD2') ,rownames(annoNew_sub)[annoNew_sub$Group_Tree=='STRd D2 Striosome']]

annoNew_sub['AB-S40304_S383_E1-50', c('marker_sum_norm_label', 'score.Tree', 'Group_Corr', 'Group_Seurat', 'Group_Tree')]

annoNew_sub[annoNew_sub$Group_Tree=='OT D1-ICj', c('marker_sum_norm_label', 'score.Tree', 'Group_Corr', 'Group_Seurat', 'Group_Tree', 'roi')]
annoNew_sub[annoNew_sub$Group_Corr=='STR SST Chodl GABA', c('marker_sum_norm_label', 'score.Tree', 'Group_Corr', 'Group_Seurat', 'Group_Tree', 'roi', 'Virus')]
annoNew_sub[annoNew_sub$Group_Corr=='STRv D1 Shell', c('marker_sum_norm_label', 'score.Tree', 'Group_Corr', 'Group_Seurat', 'Group_Tree', 'roi', 'Virus')]

length(annoNew_sub[(annoNew_sub$Group_Corr=='STRv D1 Shell')&(annoNew_sub$roi=='STRvACB'), 'Group_Tree'])
length(annoNew_sub[(annoNew_sub$Group_Corr=='STRv D1 Shell')&(annoNew_sub$roi=='STRv_ACB'), 'Group_Tree'])
table(annoNew_sub[(annoNew_sub$Group_Corr=='STRv D1 Shell')&(annoNew_sub$roi!='STRvACB')&(annoNew_sub$roi!='STRv_ACB'), 'Group_Tree'])
dim(annoNew_sub[(annoNew_sub$Group_Corr=='STRv D1 Shell')&(annoNew_sub$Group_Tree=='STRv D1 Shell'),])
table(annoNew_sub[(annoNew_sub$Group_Corr=='STRv D1 Shell')&((annoNew_sub$roi=='STRvACB')|(annoNew_sub$roi=='STRv_ACB')), 'Cluster_Corr'])
table(annoNew_sub[(annoNew_sub$Group_Corr=='STRv D1 Shell')&(annoNew_sub$roi!='STRvACB')&(annoNew_sub$roi!='STRv_ACB'), 'Cluster_Corr'])


query.data[c('SST','NPY', 'NOS1', 'TACR1', 'SCN9A'),rownames(annoNew_sub)[(annoNew_sub$Group_Corr=='STR SST Chodl GABA')&(annoNew_sub$Group_Tree=='STR SST Chodl GABA')]]
annoNew_sub[(annoNew_sub$Group_Corr=='STR SST Chodl GABA')&(annoNew_sub$Group_Tree=='STR SST Chodl GABA'), c('marker_sum_norm_label', 'score.Tree', 'Group_Corr', 'Group_Seurat', 'Group_Tree', 'roi', 'Virus')]