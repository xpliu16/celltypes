# Load doesn't work
x <- readRDS("/Users/xiaoping.liu/celltypes/NHP_BG_anal/NHP_BG_RSC_204_324_map_sub2/dend.RData")
# class(x) is dendrogram

plot(x, labels = NULL, hang=-1,
     main = "Cluster dendrogram", sub = NULL,
     xlab = NULL, ylab = "Height")

hcd = as.dendrogram(hc)
# alternative way to get a dendrogram
plot(hcd)

plot(cut(hcd, h = 0.4)$upper, main = "Upper tree")

library(ggdendro)
dend_data <- dendro_data(hcd, type = "rectangle")


library(feather)
df <- read_feather("/Users/xiaoping.liu/celltypes/NHP_BG_anal/NHP_BG_RSC_204_324_map_sub2/anno.feather")

### But the dend.RData file only has dendrograms at the cluster level, so try making it ourselves

## Get cluster medians
library(scrattch.mapping)
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114"
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping"

AIT.anndata <- read_h5ad(file.path(refFolder,"AIT_114_taxonomy.h5ad"))
# Assign levels
#subclass_lavels <- unique(AIT.anndata$uns$clusterInfo$subclass_label)

# A gene (rows) x samples (columns) sparse matrix
# A cluster factor object
cluster = AIT.anndata$uns$clusterInfo$subclass_label
names(cluster) = rownames(AIT.anndata$X)
medianExpr = get_cl_medians(t(AIT.anndata$X), cluster) ## anno should be some grouping of cells, e.g. subclass

# Get dend feature set
AIT.anndata$uns$dend = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_114/reference.rda"
load(AIT.anndata$uns$dend)
dend = reference$dend
allMarkers = unique(unlist(get_dend_markers(dend)))

##
dend.result = build_dend(
  cl.dat = medianExpr[intersect(rownames(medianExpr), allMarkers),],
  cl.cor = NULL,
  # l.color = use.color, ## Can be ignored
  # l.rank = l.rank, ## Can be ignored
  nboot = 1,
  ncores = 1)

##
dend = dend.result$dend

png(file.path(mappingFolder,'AIT_114_subclass_dend.png'), width=600)
par(mar = c(2*par("mar")[1], 1*par("mar")[2], 1*par("mar")[3], 1*par("mar")[4]))

plot(dend, labels = NULL, hang=-1,
     main = "Cluster dendrogram", sub = NULL,
     xlab = '', ylab = "Height")

dev.off()