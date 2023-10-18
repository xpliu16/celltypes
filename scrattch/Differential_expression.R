refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/"
#dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

#browseVignettes("EnhancedVolcano")

library(scrattch.mapping)
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(stringr)
library(ggplot2)
library(gridExtra)

#type1 = 'LHX6-TAC3-PLPP4'
#type2 = 'PVALB-COL19A1-ST18'
#type1 = 'MSN'
#type2 = 'IN'
# type2 = NULL

# Put this file up on server
goi = read.csv(file.path(mappingFolder,'VGIC_short.csv'))   # Genes of interest

# To see sample data from DE Seurat tutorail
#devtools::install_github('satijalab/seurat-data')
#library(SeuratData)
#InstallData("pbmc3k")
#pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")
#levels(pbmc)

#AIT.anndata <- loadTaxossh -N -L 8787:n77:8787 xiaoping.liu@hpc-loginnomy(refFolder)
# Use complete taxonomy
#AIT.anndata <- read_h5ad(file.path(refFolder, "NHP_BG_AIT115_complete.h5ad"))
AIT.anndata <- loadTaxonomy(refFolder)

annoBG <- read_feather(file.path(refFolder, "anno.feather"))

# Striatal subclasses only (at least 5% of all cells are in dSTR or vSTR)
subclasses = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
             "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
             "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
             "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
             "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
             "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   

FCcutoff = 0.5    # minimum log2 fold-change for DEG

deg_comp <- function(Expr.dat, annoBG, subclasses, goi, type1, type2, level, FCcutoff) {

  #load(paste(refFolder,"/patchseq/QC_markers.RData",sep = ""))
  #annoBG <- annoBG[match(Expr.dat$sample_id,annoBG$sample_id),]
  #annoBG <- annoBG[match(rownames(Expr.dat),annoBG$sample_id),]
  annoBG <- annoBG[match(rownames(Expr.dat),annoBG$cellNames_label),] 

  kpsubclass<-is.element(annoBG$level3.subclass_label,subclasses)
  annoBG<-annoBG[kpsubclass,]
  Expr.dat<-Expr.dat[kpsubclass,]
  #anno_type = annoBG$level3.subclass_label
  anno_type = annoBG[[level]]
  anno_type_color = annoBG[[level %>% str_replace("label", "color")]]

  # Optional: subset to ion channel genes
  genesSamp1 <- is.element(colnames(Expr.dat),goi$Approved.symbol)
  Expr.dat <- Expr.dat[,genesSamp1]
  # Missing genes
  # goi$Approved.symbol[!is.element(goi$Approved.symbol, colnames(Expr.dat)[genesSamp1])]
  # "CACNG4" "KCNF1"  "KCNK15" "KCNE1"  "CLCNKA" "KCNJ11" "KCNJ12" "KCNJ18" "ANO2"   "CNGB1" 
  # For the most part these are not showing up in Ensembl for Macaque

  Expr.dat<-t(Expr.dat)

  ident_vec <- rep("Other", length(anno_type))
  ident_vec[is.element(anno_type, type1)] = "type1"
  if (is.null(type2)){
    ident_vec[!is.element(anno_type, type1)] = "type2"
  } else {
    ident_vec[is.element(anno_type, type2)] = "type2"
  } 

  # Normalize for overall expression level of genes between cell types
  # Note that if second group is NULL, no normalization ends up being applied
  Expr.dat_norm <- Expr.dat
  t1_mean_expr = mean(Expr.dat[,ident_vec == 'type1'])
  t2_mean_expr = mean(Expr.dat[,ident_vec == 'type2'])
  expr_ratio = t1_mean_expr/t2_mean_expr
  Expr.dat_norm[,ident_vec == 'type2'] <- Expr.dat_norm[,ident_vec == 'type2'] * expr_ratio
  # check
  # t2_mean_expr_norm = mean(Expr.dat[,annoBG$level3.subclass_label == type2])

  ## For curiosity to compare p-vals, subsample data
  #keepinds = anno_type == type1 | anno_type == type2
  ##sample_factor = 10
  ##keepinds = sample(which(keepinds), round(sum(keepinds)/sample_factor))
  #Expr.dat_norm <- Expr.dat_norm[,keepinds]   # Or Expr.dat if not normalized
  #annoBG <- annoBG[keepinds,]
  #anno_type <- anno_type[keepinds]
  #anno_type_color <-anno_type_color[keepinds]

  dataBG_all<-Expr.dat_norm
  annoBG_all<-annoBG

  #brain.data     <- cbind(dataBG_all[keepGenes,],dataBG_all_PS[keepGenes,])  # Include only genes subsetted above
  brain.data <- dataBG_all
  brain.metadata <- data.frame(subclass = anno_type,
                              subclass_color = anno_type_color,
            area = annoBG_all$roi_label)
  rownames(brain.metadata) <- colnames(brain.data)

  ## Construct data set lists
  brain      <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)
  #Idents(brain) <- brain.metadata$subclass
  Idents(brain) <- ident_vec
  #brain_log <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
  de.markers <- FindMarkers(brain, ident.1 = "type1", ident.2 = "type2")

  expr1 = mean(brain.data['KCNH7', anno_type==type1])
  expr2 = mean(brain.data['KCNH7', anno_type==type2])

  log2fc_manual = log2(1+expr1)-log2(1+expr2)
  return (de.markers)
}

Expr.dat <- AIT.anndata$layers['counts']
allMarkers = c()
comparisons = list(list('MSN', NULL, 'level1.class_label', 'MSN vs rest'), 
list('PVALB-COL19A1-ST18', NULL, 'level3.subclass_label', 'FS IN vs rest'), 
list('D1-Matrix', 'D2-Matrix', 'level3.subclass_label', 'D1-Matrix vs D2-Matrix'),
list(c("D1-Matrix", "D2-Matrix"), c("D1-Striosome", "D2-Striosome"), 'level3.subclass_label', 'Matrix vs Striosome'),
list(c("D1-Matrix", "D2-Matrix"), c("D1-ShellOT", "D2-ShellOT"), 'level3.subclass_label', 'Dorsal vs Ventral'),
list('D1D2-Hybrid', c("D1-Matrix", "D2-Matrix", "D1-Striosome", "D2-Striosome"), 'level3.subclass_label', 'D1D2-Hybrid vs other dorsal MSN'),
list(c("LHX6-TAC3-PLPP4","TAC3-LHX8-PLPP4"), NULL, 'level3.subclass_label', 'TAC3-PLPP4 vs rest'))
list("SST_Chodl", NULL, 'level3.subclass_label', 'SST Chodl vs rest')
#comparisons = list(list('D1-Matrix', 'D2-Matrix', 'level3.subclass_label'))
#Expr.dat <- AIT.anndata$X   # These are already log transformed
#Expr.dat <- AIT.anndata$layers['UMIs']   # For complete taxonomy (but doesn't match anno)

for (comp in comparisons){
  print(paste0(comp[1], ' versus ', comp[2]))
  markers = deg_comp(Expr.dat, annoBG, subclasses, goi, comp[[1]], comp[[2]], comp[[3]], FCcutoff)
  print(markers)
  temp <- markers$avg_log2FC
  allMarkers <- append(allMarkers, rownames(markers)[temp>FCcutoff])
  print(allMarkers)
  allMarkers <- append(allMarkers, rownames(markers)[temp<(-FCcutoff)])
  
  #if (is.null(comp[[2]])){
  #  title = paste0(comp[[1]], collapse=' and ')
  #} else {
  #  title = paste0(paste0(comp[[1]], collapse=' and '), '_vs_', paste0(comp[[2]], collapse=' and '))
  #}
  title = comp[[4]]
  write.csv(markers, file.path(mappingFolder, paste0('DEG/', title, "_ion_channels_norm.csv")))
  png(file.path(mappingFolder, paste0('DEG/', title, "_ion_channels_norm.png")), width = 600, height = 500)
  print(EnhancedVolcano(markers,
    lab = rownames(markers),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    pCutoff = 10e-4,
    FCcutoff = FCcutoff,
  #  xlim = c(-5.5, 5.5),
  #  ylim = c(0, -log10(10e-12)),
    pointSize = 1.5,
    labSize = 2.5,
    title = title,
    subtitle = 'Ion channel DEGs, normalized',
    caption = 'FC cutoff, 0.5; p-value-adj cutoff, 10e-4',
    legendPosition = "right",
    legendLabSize = 14,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.9,
    drawConnectors = TRUE,
    typeConnectors = 'open',
    #lengthConnectors = unit(0.05, 'npc'), 
    hline = c(10e-8),
    #widthConnectors = 0.5
    )
  )
  dev.off()
}
# NOTE png's can be slow to save
allMarkers <- unique(allMarkers)

type1 = 'D1-Matrix'
type2 = 'D2-Matrix'
keepinds = anno_type == type1 | anno_type == type2
Expr.dat <- Expr.dat[,keepinds]
Expr.dat <- t(Expr.dat)
anno_type <- anno_type[keepinds]
anno_type_color <-anno_type_color[keepinds]
#genes = c("KCNC2", "KCNQ5", "KCNIP4", "CACNA2D3")
genes = allMarkers
p = list()
for (i in 1:length(genes)){
    gene = genes[i]
    #scale_fill_brewer()
    df <- data.frame(expr = Expr.dat[,gene], subclass = anno_type)
    df <- mutate(df, expr_log = log2(expr + 1)) 

    # Basic violin plot
    p[[i]] <- ggplot(df, aes(x=subclass, y=expr_log)) + geom_violin(bw = 0.3) + geom_boxplot(width=0.1) 
    p[[i]] <- p[[i]] + ggtitle(gene)
}    
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol=4)
do.call(grid.arrange,c(p, ncol=4))

# PRETTY UP COLORS
# GET DOCALL TO WORK
# COMPARE MSN vs non-MSNlibrary

#BiocManager::install("monocle")
library(monocle)

#Expr.dat <- AIT.anndata$layers['counts']
subclass_short = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome",
             "PVALB-COL19A1-ST18", "CCK-VIP-TAC3", "CCK-FBXL7", "SST_Chodl", 
             "D2-ShellOT", "CHAT", "D1-ShellOT", "D1-NUDAP", "TAC3-LHX8-PLPP4") 

Expr.sub <- Expr.dat[,allMarkers]
Expr.sub <- Expr.sub[is.element(annoBG$level3.subclass_label,subclass_short),]
Expr.sp <- t(as.matrix(Expr.sub))
Expr.sp <- as(Expr.sp, "sparseMatrix") 
fd <- data.frame(gene_short_name=rownames(Expr.sp))
annoBG <- as.data.frame(annoBG[is.element(annoBG$level3.subclass_label,subclass_short),])
rownames(annoBG) = colnames(Expr.sp) 
rownames(fd) = rownames(Expr.sp)
cds <- newCellDataSet(Expr.sp,
                         phenoData = new("AnnotatedDataFrame", data = annoBG),
                         featureData = new("AnnotatedDataFrame", data = fd))


#cds_subset <- cds[row.names(subset(Expr.sp,
#                 gene_short_name %in% allMarkers)),]
plot_genes_violin(cds, min_expr=0.1)

# group_cells_by=annoBG["level3.subclass_Tree"], 