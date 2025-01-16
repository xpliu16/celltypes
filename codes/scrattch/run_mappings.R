suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)
library (dplyr)
library(stringr)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")

run_mappings <- function(refFolder, mappingFolder, data_dir, data_fn, mode, 
                         h5ad_fn = NULL, hierarchy, proj_strs, roi_strs, 
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
    AIT.anndata$obs['off_target'] = AIT.anndata$obs[hierarchy[2]]
    
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
                                        label.cols = c(neigh_colname, class_colname, subclass_colname, cluster_colname)
  )
  a <- strsplit(refFolder,'/')[[1]]
  taxname <- a[length(a)]
  b <- strsplit(data_fn, '_')[[1]]
  dataname <- b[2]
  
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