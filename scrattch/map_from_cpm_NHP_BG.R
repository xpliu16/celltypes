#refFolder <- "/home/xiaoping.liu/scrattch/reference/NHP_BG_AIT_114"
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping"    # MAKE THIS MORE SPECIFIC
dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

suppressPackageStartupMessages({
  library("scrattch.mapping")
})
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)  # Can adjust this value if needed, depending on number of cells
options(future.rng.onMisuse="ignore")

library(scrattch.mapping)

# This remakes taxonomy from feather files and reloads, very slow, use read_h5ad instead for the time being
#AIT.anndata <- loadTaxonomy(refFolder = refFolder)

AIT.anndata <- read_h5ad(file.path(refFolder,"AIT_114_taxonomy.h5ad"))
#annoReference   = feather(file.path(refFolder,"anno.feather")) 
#exprReference   = feather(file.path(refFolder,"data.feather"))
#annoReference = as.data.frame(annoReference[match(exprReference[[sample_id]], annoReference[[sample_id]]),])
#rownames(annoReference) = rownames(datReference) = annoReference[[sample_id]]

AIT.anndata$uns$dend = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114/reference.rda"

# Was "/20220907_RSC-204-310_macaque_patchseq_star2.7_cpm.Rdata"
"/20221215_RSC-204-318_macaque_patchseq_star2.7_cpm.Rdata"
load(paste0(data_dir, "/20230309_RSC-204-324_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20230309_RSC-204-324_macaque_patchseq_star2.7_samp.dat.Rdata"))

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
load(AIT.anndata$uns$dend)
dend = reference$dend
allMarkers = unique(unlist(get_dend_markers(dend)))

## Subset query data to just those markers
query.data = query.data[intersect(rownames(query.data), allMarkers),]

# Confirm
#query.metadata$exp_component_name[1:5]
#colnames(counts)[1:5]

# Filter out unclassifieds
#kp          <- annotations$broad_type!="Unclassified"
#counts      <- counts[,kp]
#annotations <- annotations[kp,]

query.mapping <- taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.data, 
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = FALSE, 
                                  label.cols = c("class_label", "subclass_label", "cluster_label")  # Currently just use subclass_label and cluster. Columns to map against
)

#label.cols = c("class", "subclass", "supertype_label", "cluster") 

write.csv(query.mapping, file.path(mappingFolder,"NHP_BG_RSC_204_324_mapping.csv"), row.names=FALSE)
save(query.mapping, file=file.path(mappingFolder,"NHP_BG_RSC_204_324_mapping.Rdata"))

# Variable renaming
#clusters  <- unique(query.mapping$cluster)   
clusters <- unique(AIT.anndata$uns$clusterInfo$cluster)
#subclass_levels <- unique(query.mapping$subclass_Corr)
subclass_lavels <- unique(AIT.anndata$uns$clusterInfo$subclass_label)

#annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = "exp_component_name", all=TRUE) 
# Alt hack:
#rownames(query.mapping) = rownames(query.metadata)
#annotations_mapped <- merge(x = query.mapping, y = query.metadata, by.x = 0, by.y = 0, all=TRUE) 
annotations_mapped$clusters <- factor(annotations_mapped$cluster_Tree, levels=clusters)
annotations_mapped$class   <- annotations_mapped$class_Tree
annotations_mapped$class[!is.element(annotations_mapped$class,c("IN","MSN"))] = "Non-neuronal"
annotations_mapped <- annotations_mapped[match(rownames(query.metadata),annotations_mapped$exp_component_name),]   # merge resorts things
type_counts_Corr = table(annotations_mapped$subclass_Corr)
type_counts_Tree = table(annotations_mapped$subclass_Tree)

# leaf node names, each also has an ID, 0 appears to be unclassified, 
# In this taxonomy there are 49 types - max(annotations$primary_type_id)
# Filters the amount of metadata to share with taxonomy - don't worry about this for now

annotations_mapped$cluster <- factor(query.mapping$cluster, levels=clusters)  # Make into discrete levels
inds1 = grepl("STR",annotations_mapped$roi)
inds2 = annotations_mapped$library_prep_pass_fail == "Pass"
inds3 = annotations_mapped$Genes.Detected >= 1000
inds4 = annotations_mapped$percent_reads_aligned_total >= 50
anno_mapped_sub = annotations_mapped[inds1&inds2&inds3&inds4,]
query.data_sub = query.data[,inds1&inds2&inds3&inds4]

type_counts_Corr_QC = table(anno_mapped_sub$subclass_Corr)
type_counts_Tree_QC = table(anno_mapped_sub$subclass_Tree)

write.csv(annotations_mapped, file.path(mappingFolder,"NHP_BG_RSC_204_324_ann_map.csv"), row.names=FALSE)
save(annotations_mapped, anno_mapped_sub, type_counts_Corr, type_counts_Tree, type_counts_Corr_QC,
     type_counts_Tree_QC, file=file.path(mappingFolder,"NHP_BG_RSC_204_324_ann_map.Rdata"))

D1_expr <- query.data_sub["DRD1",]
D1_expr_D1_types = D1_expr[is.element(anno_mapped_sub$subclass_Tree,c("D1-Matrix","D1-ShellOT","D1-Striosome", "D1-ICj"))]
D1_expr_D1_types <- na.omit(D1_expr_D1_types)
length(D1_expr_D1_types)
D1_expr_D1_types = mean(D1_expr_D1_types)
D1_expr_D2_types = D1_expr[is.element(anno_mapped_sub$subclass_Tree,c("D2-Matrix","D2-ShellOT","D2-Matrix / D2-Striosome", "D2-Striosome"))]
D1_expr_D2_types <- na.omit(D1_expr_D2_types)
D1_expr_D2_types = mean(D1_expr_D2_types)

D2_expr <- query.data_sub["DRD2",]
D2_expr_D1_types = D2_expr[is.element(anno_mapped_sub$subclass_Tree,c("D1-Matrix","D1-ShellOT","D1-Striosome", "D1-ICj"))]
D2_expr_D1_types <- na.omit(D2_expr_D1_types)
length(D2_expr_D1_types)
D2_expr_D1_types = mean(D2_expr_D1_types)
D2_expr_D2_types = D2_expr[is.element(anno_mapped_sub$subclass_Tree,c("D2-Matrix","D2-ShellOT","D2-Matrix / D2-Striosome", "D2-Striosome"))]
D2_expr_D2_types <- na.omit(D2_expr_D2_types)
D2_expr_D2_types = mean(D2_expr_D2_types)

RXFP1_expr <- query.data_sub["RXFP1",]
RXFP1_expr_D1D2_types = RXFP1_expr[is.element(anno_mapped_sub$subclass_Tree,c("D1D2 Hybrid"))]
#RXFP1_expr_D1D2_types <- na.omit(RXFP1_expr_D1D2_types)
length(RXFP1_expr_D1D2_types)
RXFP1_expr_D1D2_types = mean(RXFP1_expr_D1D2_types)
# May exclude "D2-Matrix / D2-Striosome" which may express RXFP1, but doesn't make a difference
indsA = is.element(anno_mapped_sub$subclass_Tree, c("D1-Matrix","D1-ShellOT","D1-Striosome", "D2-Matrix / D2-Striosome", "D1-ICj", "D2-Matrix","D2-ShellOT", "D2-Striosome"))
indsB = !is.element(anno_mapped_sub$cluster_Tree, c("53_MSN", "54_MSN", "55_MSN"))
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

table(anno_mapped_sub$subclass_Tree)

#counts = counts[intersect(rownames(counts), allMarkers),]
query.counts_sub <- counts[,inds1&inds2&inds3&inds4]
query.metadata_sub <- query.metadata[inds1&inds2&inds3&inds4,]
query.mapping_sub <- query.mapping[inds1&inds2&inds3&inds4,]

dir.create(file.path(mappingFolder, "NHP_BG_RSC_204_324_map_sub2"))

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = file.path(mappingFolder, "NHP_BG_RSC_204_324_map_sub2"),
                      query.data     = query.counts_sub,  # Don't need log-normalized data here
                      query.metadata = query.metadata_sub,
                      query.mapping  = query.mapping_sub,
                      doPatchseqQC   = FALSE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

dir.create(file.path(mappingFolder, "NHP_BG_RSC_204_324_map"))

buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = file.path(mappingFolder, "NHP_BG_RSC_204_324_map"),
                      query.data     = counts,  # Don't need log-normalized data here
                      query.metadata = query.metadata,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = FALSE,  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
)

dir(mappingFolder)


anno_hybrid_tree = anno_mapped_sub[anno_mapped_sub$subclass_Tree=='D1D2 Hybrid',]


library(feather)
df <- read_feather(file.path(mappingFolder,"anno.feather"))
