#refFolder <- "/home/xiaoping.liu/scrattch/reference/tasic_2016"
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping"
dir.create(mappingFolder, showWarnings=FALSE)

suppressPackageStartupMessages({
  library("scrattch.mapping")
  library("tasic2016data")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)
library(tasic2016data)

AIT.anndata <- loadTaxonomy(refFolder = refFolder)
#load(file.path(refFolder,"reference.rda"))
#AIT.anndata <- reference

annotations <- tasic_2016_anno
counts      <- tasic_2016_counts

# Put annotations and counts in the same order
annotations <- annotations[match(colnames(counts),annotations$sample_name),] 
rownames(annotations) <- annotations$sample_name  

# Filter out unclassifieds
kp       <- subsampleCells(annotations$cluster,subSamp=5,seed=10)&(annotations$class!="Non-neuronal")
#kp          <- annotations$broad_type!="Unclassified"
counts      <- counts[,kp]
annotations <- annotations[kp,]

# Variable renaming
clusters            <- annotations$primary_type_label[match(1:49,annotations$primary_type_id)]   
# leaf node names, each also has an ID, 0 appears to be unclassified, 
# In this taxonomy there are 49 types - max(annotations$primary_type_id)
# Filters the amount of metadata to share with taxonomy - don't worry about this for now
annotations$cluster <- factor(annotations$primary_type_label, levels=clusters)   # Make into discrete levels
annotations$class   <- annotations$broad_type
# Change annotation of anything that is not GABA-ergic or Glutamatergic to Non-neuronal
annotations$class[!is.element(annotations$class,c("GABA-ergic Neuron","Glutamatergic Neuron"))] = "Non-neuronal"

query.counts   <- counts
query.metadata <- annotations
query.logCPM   <- logCPM(query.counts)

query.mapping <- taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.logCPM, 
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = TRUE, 
                                  label.cols = c("cluster_label") # Columns to map against
)

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = mappingFolder,
                      query.data     = query.counts,  # Don't need log-normalized data here
                      query.metadata = query.metadata,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = FALSE  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

dir(mappingFolder)


