refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/"
#dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

library(scrattch.mapping)
library(Seurat)
library(dplyr)

type1 = 'LHX6-TAC3-PLPP4'
type2 = 'PVALB-COL19A1-ST18'

# Put this file up on server
goi = read.csv(file.path(mappingFolder,'VGIC_short.csv'))   # Genes of interest

# To see sample data from DE Seurat tutorail
#devtools::install_github('satijalab/seurat-data')
#library(SeuratData)
#InstallData("pbmc3k")
#pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")
#levels(pbmc)


#AIT.anndata <- loadTaxonomy(refFolder)
# Use complete taxonomy
AIT.anndata <- read_h5ad(file.path(refFolder, "NHP_BG_AIT115_complete.h5ad"))
#Expr.dat <- AIT.anndata$X
Expr.dat <- AIT.anndata$layers['UMIs']
annoBG <- read_feather(file.path(refFolder, "anno.feather"))
#load(paste(refFolder,"/patchseq/QC_markers.RData",sep = ""))
#annoBG <- annoBG[match(Expr.dat$sample_id,annoBG$sample_id),]
annoBG <- annoBG[match(rownames(Expr.dat),annoBG$sample_id),]

# Striatal subclasses only (at least 5% of all cells are in dSTR or vSTR)
subclass = c("D1-Matrix", "D2-Striosome", "D2-Matrix", "D2-Hybrid-MCHR2", 
             "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", "SLC17A7-SATB2",
             "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", "CCK-FBXL7",
             "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
             "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
             "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
             "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   


kpsubclass<-is.element(annoBG$level3.subclass_label,subclass)
annoBG<-annoBG[kpsubclass,]
Expr.dat<-Expr.dat[kpsubclass,]

# Optional: subset to ion channel genes
genesSamp1 <- is.element(colnames(Expr.dat),goi$Approved.symbol)
Expr.dat <- Expr.dat[,genesSamp1]

Expr.dat<-t(Expr.dat)

# Normalize for overall expression level of genes between cell types
Expr.dat_norm <- Expr.dat
t1_mean_expr = mean(Expr.dat[,annoBG$level3.subclass_label == type1])
t2_mean_expr = mean(Expr.dat[,annoBG$level3.subclass_label == type2])
expr_ratio = t1_mean_expr/t2_mean_expr
Expr.dat_norm[,annoBG$level3.subclass_label == type2] <- Expr.dat_norm[,annoBG$level3.subclass_label == type2] * expr_ratio
# check
# t2_mean_expr_norm = mean(Expr.dat[,annoBG$level3.subclass_label == type2])

keepinds = annoBG$level3.subclass_label == type1 | annoBG$level3.subclass_label == type2
# For curiosity to compare p-vals, subsample data
#sample_factor = 10
#keepinds = sample(which(keepinds), round(sum(keepinds)/sample_factor))
Expr.dat_norm <- Expr.dat_norm[,keepinds]   # Or Expr.dat if not normalized
annoBG <- annoBG[keepinds,]

dataBG_all<-Expr.dat_norm
annoBG_all<-annoBG

#brain.data     <- cbind(dataBG_all[keepGenes,],dataBG_all_PS[keepGenes,])  # Include only genes subsetted above
brain.data <- dataBG_all
brain.metadata <- data.frame(subclass = annoBG_all$level3.subclass_label,
                             subclass_color = annoBG_all$level3.subclass_color,
			     area = annoBG_all$roi_label)
rownames(brain.metadata) <- colnames(brain.data)

## Construct data set lists
brain      <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)
Idents(brain) <- brain.metadata$subclass
#brain_log <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
de.markers <- FindMarkers(brain, ident.1 = type1, ident.2 = type2)

write.csv(de.markers, file.path(mappingFolder, paste0("/DEG/",type1,"-vs-", type2, "_ion_channels_norm.csv")))

expr1 = mean(brain.data['KCNH7', annoBG_all$level3.subclass_label==type1])
expr2 = mean(brain.data['KCNH7', annoBG_all$level3.subclass_label==type2])

log2fc_manual = log2(1+expr1)-log2(1+expr2)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

browseVignettes("EnhancedVolcano")

library(EnhancedVolcano)

EnhancedVolcano(de.markers,
  lab = rownames(de.markers),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 10e-4,
  FCcutoff = 0.5,
#  xlim = c(-5.5, 5.5),
#  ylim = c(0, -log10(10e-12)),
  pointSize = 1.5,
  labSize = 2.5,
  title = paste0(type1, ' vs ', type2, ' ion channel genes, normalized'),
#  title = paste0(type1, ' vs ', type2, 'all genes'),
  subtitle = 'Differential expression',
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

library(ggplot2)
Expr.dat <- Expr.dat[,keepinds]
Expr.dat <- t(Expr.dat)
genes = c("KCNC2", "KCNQ5", "KCNIP4", "CACNA2D3")
for (i in 1:length(genes)){
    gene = genes[i]
    #scale_fill_brewer()
    df <- data.frame(expr = Expr.dat[,gene], subclass = annoBG_all$level3.subclass_label)
    df <- mutate(df, expr_log = log2(expr + 1)) 

    # Basic violin plot
    p[[i]] <- ggplot(df, aes(x=subclass, y=expr_log)) + geom_violin(bw = 0.3) + geom_boxplot(width=0.1) 
    p[[i]] <- p[[i]] + ggtitle(gene)
}    
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol=4)
do.call(grid.arrange,c(p, ncol=4))

# PRETTY UP COLORS
# GET DOCALL TO WORK
# COMPARE MSN vs non-MSN
