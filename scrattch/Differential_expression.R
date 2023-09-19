#refFolder <- "/home/xiaoping.liu/scrattch/reference/NHP_BG_AIT_114"
refFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115"
#mappingFolder <- paste0(refFolder,"/mapping/")
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/"
#dir.create(mappingFolder, showWarnings=FALSE)
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object"

library(scrattch.mapping)
library(Seurat)

# Put this file up on server
goi = read.csv(file.path(mappingFolder,'HUGO_genes/VGIC_short.csv'))  # Genes of interest

load(paste0(data_dir, "/20230807_RSC-204-337_macaque_patchseq_star2.7_cpm.Rdata"))
load(paste0(data_dir, "/20230807_RSC-204-337_macaque_patchseq_star2.7_samp.dat.Rdata"))

AIT.anndata <- loadTaxonomy(refFolder)
Expr.dat <- AIT.anndata$X
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

# subset to patch dend genes
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
dend <- readRDS(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
dendMarkers = unique(unlist(get_dend_markers(dend)))
genesSamp1 <- is.element(colnames(Expr.dat),dendMarkers)
Expr.dat <- Expr.dat[,genesSamp1]
genesSamp2 <- (is.element(colnames(patch.dat),dendMarkers))
patch.dat <- patch.dat[,genesSamp2]

Expr.dat<-t(Expr.dat)

dataBG_all<-Expr.dat
annoBG_all<-annoBG
dataBG_all_PS<-patch.dat
annoBG_all_PS<-patch_anno

brain.data     <- cbind(dataBG_all[keepGenes,],dataBG_all_PS[keepGenes,])  # Include only genes subsetted above
brain.metadata <- data.frame(set=c(rep("FACs",dim(dataBG_all)[2]),rep("PatchSeq",dim(dataBG_all_PS)[2])),
                             subclass = c(annoBG_all$level3.subclass_label,annoBG_all_PS$level3.subclass_Tree),
                             QC = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$Norm_Marker_Sum.0.4_label), 
                             QC2 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$NMS_stringent),
                             QC3 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$library_prep_pass_fail),
                             QC4 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$percent_reads_aligned_to_introns>25),
                             QC5 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$score.Corr>0.55),
                             QC6 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$Genes.Detected<8000),
                             QC6 = c(rep("none",dim(dataBG_all)[2]), annoBG_all_PS$amplified_quantity_ng),
                             subclass_color = c(annoBG_all$level3.subclass_color, annoBG_all$level3.subclass_color[match(annoBG_all_PS$level3.subclass_Tree, annoBG_all$level3.subclass_label)]),
                             area =c(rep("BG",dim(dataBG_all)[2]),annoBG_all_PS$Region),
                             subclass_corr =c(annoBG_all$level3.subclass_label,annoBG_all_PS$level3.subclass_Corr),
                             subclass_color_corr = c(annoBG_all$level3.subclass_color, annoBG_all$level3.subclass_color[match(annoBG_all_PS$level3.subclass_Corr, annoBG_all$level3.subclass_label)]),
                             prep = c(rep("acute",dim(dataBG_all)[2]),annoBG_all_PS$Prep))
rownames(brain.metadata) <- colnames(brain.data)

## Construct data set lists
brain      <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)

