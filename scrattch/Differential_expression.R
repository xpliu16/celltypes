#refFolderList <- list("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_115")
refFolderList <- list("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_MTG",
"/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_V1", 
"/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_M1")
#refFolderList <- list("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_V1") # For within region cluster comparison
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_115/"

#BiocManager::install("EnhancedVolcano")

#browseVignettes("EnhancedVolcano")

library(scrattch.mapping)
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(stringr)
library(ggplot2)
library(gridExtra)
library(plyr)
library(RColorBrewer)

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

# Use complete taxonomy
#AIT.anndata <- read_h5ad(file.path(refFolder, "NHP_BG_AIT115_complete.h5ad"))
#AIT.anndata<-loadTaxonomy(refFolder)
anndataList <- list()    # Raw data
#annoList <- list()         # Metadata

for (refFolder in refFolderList) {
  pathspl <- strsplit(refFolder,'/')[[1]]
  annName <- pathspl[length(pathspl)]
  anndataList[[annName]] <- loadTaxonomy(refFolder)
  #annoList[[annName]] <- read_feather(file.path(refFolder, "anno.feather"))
}

FCcutoff = 0.5    # minimum log2 fold-change for DEG

deg_comp <- function(Expr.dat, anno, goi, type1, type2, level, FCcutoff, 
                    autocolor = FALSE, normalize = FALSE) {

  #load(paste(refFolder,"/patchseq/QC_markers.RData",sep = ""))
  #anno <- anno[match(Expr.dat$sample_id,anno$sample_id),]
  anno <- anno[match(rownames(Expr.dat),anno$sample_id),]
  #anno <- anno[match(rownames(Expr.dat),anno$cellNames_label),] 

  #kpsubclass<-is.element(anno$level3.subclass_label,subclasses)
  #anno<-anno[kpsubclass,]
  #Expr.dat<-Expr.dat[kpsubclass,]
  #anno_type = anno$level3.subclass_label
  print(level)
  anno_type = anno[[level]]
  print(anno_type)
  if (autocolor) {
    anno_type_color = anno[[level %>% str_replace("label", "color")]]
  } else {
    types <- unique(anno_type)
    pal <- brewer.pal(n = length(types), name = "Set2")[1:length(types)] # Min length is 3, so need to crop
    anno_type_color <- mapvalues(anno_type, types, pal)
  }
  # Optional: subset to ion channel genes
  genesSamp1 <- is.element(colnames(Expr.dat),goi$Approved.symbol)
  Expr.dat <- Expr.dat[,genesSamp1]
  # Missing genes
  # goi$Approved.symbol[!is.element(goi$Approved.symbol, colnames(Expr.dat)[genesSamp1])]
  # "CACNG4" "KCNF1"  "KCNK15" "KCNE1"  "CLCNKA" "KCNJ11" "KCNJ12" "KCNJ18" "ANO2"   "CNGB1" 
  # For the most part these are not showing up in Ensembl for Macaque

  Expr.dat<-t(Expr.dat)

  print(type1)
  print(type2)
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
  if (normalize) {
    t1_mean_expr = mean(Expr.dat[,ident_vec == 'type1'])
    t2_mean_expr = mean(Expr.dat[,ident_vec == 'type2'])
    expr_ratio = t1_mean_expr/t2_mean_expr
    Expr.dat_norm[,ident_vec == 'type2'] <- Expr.dat_norm[,ident_vec == 'type2'] * expr_ratio
    # check
    # t2_mean_expr_norm = mean(Expr.dat[,anno$level3.subclass_label == type2])
  } 

  ## For curiosity to compare p-vals, subsample data
  #keepinds = anno_type == type1 | anno_type == type2
  ##sample_factor = 10
  ##keepinds = sample(which(keepinds), round(sum(keepinds)/sample_factor))
  #Expr.dat_norm <- Expr.dat_norm[,keepinds]   # Or Expr.dat if not normalized
  #anno <- anno[keepinds,]
  #anno_type <- anno_type[keepinds]
  #anno_type_color <-anno_type_color[keepinds]

  dataBG_all<-Expr.dat_norm
  anno_all<-anno

  #brain.data     <- cbind(dataBG_all[keepGenes,],dataBG_all_PS[keepGenes,])  # Include only genes subsetted above
  brain.data <- dataBG_all
  brain.metadata <- data.frame(type = anno_type,
                              type_color = anno_type_color)
  #                            area = anno_all$roi_label)
  rownames(brain.metadata) <- colnames(brain.data)

  ## Construct data set lists
  brain      <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)
  #Idents(brain) <- brain.metadata$subclass
  Idents(brain) <- ident_vec
  #brain_log <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
  if (length(unique(ident_vec))==1) {
    return (NULL)
  } else {
    de.markers <- FindMarkers(brain, ident.1 = "type1", ident.2 = "type2")    # If you want to downsample: max.cells.per.ident = max_n_cells
  
    expr1 = mean(brain.data['KCNH7', anno_type==type1])
    expr2 = mean(brain.data['KCNH7', anno_type==type2])

    log2fc_manual = log2(1+expr1)-log2(1+expr2)
    return (de.markers)
  }
}

allMarkers <- c()
Expr.dat.List <- list()
annoAll <- data.frame()
#subset=NULL   # ~48640 is ngenes, slightly different between taxonomies
# Striatal subclasses only (at least 5% of all cells are in dSTR or vSTR)
#subset = list('level3.subclass_label', c("D1-Matrix", "D2-Striosome", "D2-Matrix", 
#             "D2-Hybrid-MCHR2", "D1D2-Hybrid", "LHX6-TAC3-PLPP4", "D1-Striosome", 
#             "SLC17A7-SATB2", "PVALB-COL19A1-ST18", "LHX6-SATB1", "CCK-VIP-TAC3", 
#             "CCK-FBXL7", "SST-RSPO2", "SST_Chodl", "D2-ShellOT", "CHAT", "D1-ShellOT",
#             "D1-NUDAP", "TAC3-LHX8-PLPP4", "MEIS2", "SN_STH_GPe-MEIS2-OTX2",
#             "LHX6-LHX8-GBX1", "LHX6_SST", "NAc-CCK-SEMA3A", "GP-LHX6",
#             "SST-ADARB2", "SLC17A6", "WDR49-ADAM12", "D1-ICj", "NAc-LHX8")   

#subset = list('CrossArea_cluster_label', c('L2/3 IT_1', 'L2/3 IT_2', 'L2/3 IT_3', 'L2/3 IT_4', 
#                                           'L2/3 IT_5', 'L2/3 IT_6', 'L2/3 IT_7','L2/3 IT_8',
#                                           'L2/3 IT_9','L2/3 IT_10','L2/3 IT_11','L2/3 IT_12','L2/3 IT_13'))

#subset = list('CrossArea_subclass_label', c('L2/3 IT'))
#subset = list('CrossArea_subclass_label', c('L4 IT'))
subset = list('CrossArea_subclass_label', c('L2/3 IT', 'L4 IT'))
# subclass = L2/3 IT includes L4, L5, and L6 IT clusters
genes = list()
for (ann in names(anndataList)) {
  AIT.anndata = anndataList[[ann]]
  Expr.dat <- AIT.anndata$layers['counts']
  anno = AIT.anndata$obs
  if (!is.null(subset)){
    subInds = is.element(anno[[subset[[1]]]], subset[[2]])
    Expr.dat <- Expr.dat[subInds,]
    anno <- anno[subInds,]  
  }
  print(dim(Expr.dat))
  print(dim(anno))
  Expr.dat.List[[ann]] <- Expr.dat
  anno$ann_source = ann
  annoAll <- rbind(annoAll, anno)
  genes[[ann]] = colnames(Expr.dat)
}
genes = reduce(genes, intersect)
Expr.dat.df <- data.frame()
for (ann in names(anndataList)) {
  Expr.dat.df <- rbind(Expr.dat.df, Expr.dat.List[[ann]][,genes])
} 

comparisons = list(list(#type1 = 'CrossAreal_V1', 
                        type1 = names(anndataList)[1],
                        type2 = names(anndataList)[2],
                        colname = 'ann_source', 
                        title = 'L23IT:_V1_vs_MTG',
                        splitby = 'CrossArea_cluster_label')) 

comparisons = list(list(type1 = 'L2/3 IT_1', 
                        type2 = c('L2/3 IT_2','L2/3 IT_5','L2/3 IT_6'), 
                        colname = 'CrossArea_cluster_label', 
                        title = 'L23IT_1_vs_L23IT_2_5_6',
                        splitby = 'ann_source'),
                   list(type1 = 'L2/3 IT_5', 
                   type2 = c('L2/3 IT_1','L2/3 IT_2','L2/3 IT_6'),
                   colname = 'CrossArea_cluster_label', 
                   title = 'L23IT_5_vs_L23IT_1_2_6',
                   splitby = 'ann_source'))

#comparisons = list(list('MSN', NULL, 'level1.class_label', 'MSN vs rest'), 
#list('PVALB-COL19A1-ST18', NULL, 'level3.subclass_label', 'FS IN vs rest'), 
#list('D1-Matrix', 'D2-Matrix', 'level3.subclass_label', 'D1-Matrix vs D2-Matrix'),
#list(c("D1-Matrix", "D2-Matrix"), c("D1-Striosome", "D2-Striosome"), 'level3.subclass_label', 'Matrix vs Striosome'),
#list(c("D1-Matrix", "D2-Matrix"), c("D1-ShellOT", "D2-ShellOT"), 'level3.subclass_label', 'Dorsal vs Ventral'),
#list('D1D2-Hybrid', c("D1-Matrix", "D2-Matrix", "D1-Striosome", "D2-Striosome"), 'level3.subclass_label', 'D1D2-Hybrid vs other dorsal MSN'),
#list(c("LHX6-TAC3-PLPP4","TAC3-LHX8-PLPP4"), NULL, 'level3.subclass_label', 'TAC3-PLPP4 vs rest'),
#list("SST_Chodl", NULL, 'level3.subclass_label', 'SST Chodl vs rest'),
#list("CHAT", NULL, 'level3.subclass_label', 'CHAT vs rest'))
#comparisons = list(list('D1-Matrix', 'D2-Matrix', 'level3.subclass_label'))
#Expr.dat <- AIT.anndata$X   # These are already log transformed
#Expr.dat <- AIT.anndata$layers['UMIs']   # For complete taxonomy (but doesn't match anno)
autocolor = FALSE
normalize = TRUE

for (comp in comparisons){
  if (!is.null(comp$splitby)){
    spl_types = unique(as.vector(annoAll[[comp$splitby]]))
    print(spl_types)
  } else {
    spl_types = 'notsplit'
  }
  
  for (subtype in spl_types) {
    if (!(subtype=='notsplit')) {
      # Subset taxonomy further
      print(subtype)
      annoAll.sub = annoAll[annoAll[comp$splitby] == subtype,]
      Expr.dat.sub = Expr.dat.df[annoAll[comp$splitby] == subtype,]
    } else {
      Expr.dat.sub <- Expr.dat.df
      annoAll.sub <- annoAll
      print('here')
    }
    print(paste0(comp$type1, ' versus ', comp$type2))
    markers = deg_comp(Expr.dat.sub, annoAll.sub, goi, comp$type1, comp$type2, 
                       comp$colname, FCcutoff, autocolor, normalize)
    if (!is.null(markers)) {
      print(markers)
      temp <- markers$avg_log2FC
      allMarkers <- append(allMarkers, rownames(markers)[temp>FCcutoff])
      print(allMarkers)
      allMarkers <- append(allMarkers, rownames(markers)[temp<(-FCcutoff)])
    
      title = paste(comp$title, gsub('/| ','',subtype), sep='_')
      if (normalize) {
        title = paste(title, 'norm', sep='_')
      }
      if (title == 'CHAT vs rest') {
        caption = 'FC cutoff, 0.5; p-value-adj cutoff, 10e-4'
      } else {
        caption = ''
      }
      write.csv(markers, file.path(mappingFolder, paste0('DEG/', title, "_ion_channels.csv")))
      png(file.path(mappingFolder, paste0('DEG/', title, "_ion_channels.png")), width = 800, height = 600)
      p1 <- EnhancedVolcano(markers,
        lab = rownames(markers),
        x = 'avg_log2FC',
        y = 'p_val_adj',
        #pCutoff = 10e-4,
        pCutoff = 0.01,
        FCcutoff = FCcutoff,
      #  xlim = c(-5.5, 5.5),
      #  ylim = c(0, -log10(10e-12)),
        pointSize = 1.5,
        labSize = 5,
        title = title,
        #subtitle = 'Ion channel DEGs, normalized',
        subtitle = 'Ion channel DEGs',
        caption = caption,
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
      xr = layer_scales(p1)$x$get_limits() %>% range
      xlim(xr*0.9)
  
      p1 <- p1 + theme(axis.text.x = element_text(size=26), axis.text.y=element_text(size=26), axis.title=element_text(size=26))
    
      print(p1)
      dev.off()
    }
  }
}
# NOTE png's can be slow to save
allMarkers <- unique(allMarkers)

library(lme4)
library (lmerTest) 
# From lmerTest: lmer overloads lme4::lmer and produced an object of class lmerModLmerTest which inherits
# from lmerMod. In addition to computing the model (using lme4::lmer), lmerTest::lmer
# computes a couple of components needed for the evaluation of Satterthwaite’s denominator
# degrees of freedom.
#gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#             data = cbpp, family = binomial))

# Normalization - this should happen before subsetting genes
rn <- rownames(Expr.dat.df)
cn <- colnames(Expr.dat.df)
Expr.dat.cpm <- as.data.frame(t(cpm(t(Expr.dat.df))))
# CPM takes (or coerces to) a matrix 
# with features / genes as rows, samples / cells as columns
rownames(Expr.dat.cpm) <- rn
colnames(Expr.dat.cpm) <- cn
# Poisson wants counts - gets rid of warnings but very small effect on regression output
Expr.dat.cpm <- round(Expr.dat.cpm)  

# Subset to ion channel genes
Expr.dat.goi <- Expr.dat.cpm[,goi$Approved.symbol]

# Remove genes expressed globally at low levels
# Expression > 2 CPM in at least min_samps, 
# where min_samps is half the size of smallest cluster
temp  = table(annoAll$CrossArea_cluster_label)
temp = temp[temp !=0]
min_samps = round(0.5 * min(temp))
isexpr <- colSums(Expr.dat.goi > 2) >= min_samps
Expr.dat.goi <- Expr.dat.goi[,isexpr]
Expr.dat.log <- log2(Expr.dat.goi +1)

#subclass = 'L2/3 IT'
subclass = 'L4 IT'

# Make sure these are the same order
all(rownames(Expr.dat.goi) == rownames(annoAll))
#all(rownames(Expr.dat.df) %in% rownames(annoAll)

#gene_set = c('KCNQ5','DPP10')
#gene_set = c('KCNQ5', 'DPP10', 'FILIP1')
#model_data = cbind(Expr.dat.goi[gene_set],
#model_data = cbind(Expr.dat.goi,
model_data = cbind(Expr.dat.log,
                   annoAll['ann_source'], 
                   annoAll['CrossArea_cluster_label'], 
                   as.factor(annoAll$Donor_id))   # This was not a factor although numeric
colnames(model_data) = c(colnames(Expr.dat.goi), 'region', 'cluster', 'donor')

# But you have no random effects
#gm1 <- glmer(kcnq5 ~ region + cluster, family = poisson, data = model_data,
#             subset = annoAll['CrossArea_subclass_label'] == subclass)
# data: an optional data frame containing the variables named in formula
# subset:	an optional expression indicating the subset of the rows of data that 
# should be used in the fit. This can be a logical vector, or a numeric vector 
# indicating which observation numbers are to be included, or a character 
# vector of the row names to be included.

#mod<-glm(KCNQ5 ~ region + cluster + region:cluster, family = poisson, data = model_data,
#                 subset = annoAll['CrossArea_subclass_label'] == subclass)
# No additional subsetting
#mod<-glm(KCNQ5 ~ region + cluster + region:cluster, family = poisson, data = model_data)
# Expands to same as above
#mod<-glm(KCNQ5 ~ region*cluster, family = poisson, data = model_data,
#         subset = annoAll['CrossArea_subclass_label'] == subclass)
coeffs = list()
pvals = list()
regional = list()  # the whole L2/3 IT and L4 IT group shows an overall effect of region on expression 
cluster = list()   # the gene is differential across clusters
reg_clust = list() # the gene is differential for a cluster that is highly regional in its expression pattern
regxclust = list() # the gene shows an additional effect of region on the specific cluster
chisq = list()
chisq_regional = list() # gene shows significant difference of adding region as a regressor
chisq_cluster = list() # gene shows significant difference of adding cluster as a regressor 
Pr_F = list()
Pr_regional = list()
Pr_cluster = list()
donor_ef = list()

psig = 0.01
coeff_thresh = 0.5
#coeff_thresh = 0.4  # for glmer

#for (gene in c('KCNQ5','DPP10')){
for (gene in colnames(Expr.dat.log)){
  gene <- as.symbol(gene)
  #mod<-eval(bquote(glm(.(gene) ~ region*cluster, family = poisson, data = model_data)))
  #mod<-eval(bquote(glmer(.(gene) ~ region*cluster + (1|donor), family = poisson, data = model_data, nAGQ = 0))) # These pvalues are Wald tests of the model with and without each individual parameter (but not that whole categorical variable)
  mod<-eval(bquote(lmer(.(gene) ~ region*cluster + (1|donor), data = model_data)))
  #mod<-eval(bquote(lme(.(gene) ~ region*cluster, random=~1|donor, data = model_data)))  # this can't handle rank deficiency from not all combinations of crossed variable having samples
  summary(mod)  
  #mod$coefficients
  #mod$residuals
  coeffs <- c(coeffs, list(coef(summary(mod))[,1]))  # coefficient estimates
  pval_i<-coef(summary(mod))[,5]  # column 4 for glm()
  pvals <- c(pvals, list(pval_i))  # pvals
  
  regional <- c(regional, 
                any((abs(coef(summary(mod))[c('regionCrossAreal_MTG',
                                             'regionCrossAreal_V1'),1]) > coeff_thresh)
                    & (p.adjust(pval_i[c('regionCrossAreal_MTG',
                                           'regionCrossAreal_V1')], 
                                method = "bonferroni") < psig)))
  
  cluster <- c(cluster, 
                any((abs(coef(summary(mod))[c('clusterL2/3 IT_5', 
                                             'clusterL2/3 IT_6',
                                             'clusterL4 IT_1',
                                             'clusterL4 IT_5',
                                             'clusterL2/3 IT_2',   
                                             'clusterL4 IT_2', 
                                             'clusterL4 IT_3', 
                                             'clusterL4 IT_4',
                                             'clusterL4 IT_6',
                                             'clusterL2/3 IT_3', 
                                             'clusterL2/3 IT_4'),1]) > coeff_thresh)
                    & (p.adjust(pval_i[c('clusterL2/3 IT_5', 
                                            'clusterL2/3 IT_6',
                                            'clusterL4 IT_1',
                                            'clusterL4 IT_5',
                                            'clusterL2/3 IT_2',   
                                            'clusterL4 IT_2', 
                                            'clusterL4 IT_3', 
                                            'clusterL4 IT_4',
                                            'clusterL4 IT_6',
                                            'clusterL2/3 IT_3', 
                                            'clusterL2/3 IT_4')],
                                method = "bonferroni") < psig)))
  reg_clust <- c(reg_clust, 
                 any((abs(coef(summary(mod))[c('clusterL2/3 IT_2', 
                                               'clusterL2/3 IT_3',
                                               'clusterL2/3 IT_4',
                                               'clusterL4 IT_2',
                                               'clusterL4 IT_3', 
                                               'clusterL4 IT_4',
                                               'clusterL4 IT_6'),1]) > coeff_thresh)
                     & (p.adjust(pval_i[c('clusterL2/3 IT_2', 
                                             'clusterL2/3 IT_3',
                                             'clusterL2/3 IT_4',
                                             'clusterL4 IT_2',
                                             'clusterL4 IT_3', 
                                             'clusterL4 IT_4',
                                             'clusterL4 IT_6')],
                                 method = "bonferroni") < psig)))
  
  terms = rownames(summary(mod)$coefficients)
  terms_sub = c(terms[grepl("regionCrossAreal_MTG:cluster", terms)],terms[grepl("regionCrossAreal_V1:cluster", terms)])
  regxclust <- c(regxclust,
                 any((abs(coef(summary(mod))[c(terms_sub),1]) > 0.4)
                     & (p.adjust(pval_i[c(terms_sub)], method = "bonferroni") < psig)))
  
  #for glm():
  ##anova(mod) # These are sequentially added 
  #chisq_i<- anova(mod,test='Chisq')['Pr(>Chi)']
  #chisq<- c(chisq, list(chisq_i))
  #chisq_regional<- c(chisq_regional, p.adjust(chisq_i['region', 'Pr(>Chi)'], method = "bonferroni") < psig)
  #chisq_cluster<- c(chisq_cluster, p.adjust(chisq_i['cluster', 'Pr(>Chi)'], method = "bonferroni") < psig)
  
  Pr_i <- anova(mod)['Pr(>F)']
  Pr_F <- c(Pr_F, list(Pr_i))
  Pr_regional <- c(Pr_regional, p.adjust(Pr_i['region', 'Pr(>F)'], method = "bonferroni") < psig)
  Pr_cluster <- c(Pr_cluster, p.adjust(Pr_i['cluster', 'Pr(>F)'], method = "bonferroni") < psig)
  donor_ef <- c(donor_ef, ranef(mod)$donor)
}  

mod_df <- data.frame(row.names = colnames(Expr.dat.log), 
                     regional = unlist(regional), 
                     cluster = unlist(cluster), 
                     reg_cluster = unlist(reg_clust), 
                     regxclust = unlist(regxclust),
                     Pr_regional = unlist(Pr_regional),
                     Pr_cluster = unlist(Pr_cluster))
                     #chisq_regional = unlist(chisq_regional),
                     #chisq_cluster = unlist(chisq_cluster))
mod_df <- cbind(glm_df,t(do.call(cbind, Pr_F)), do.call(rbind,pvals))
# If probabilities are not significant, blank out the coefficient as NA
temp <- do.call(rbind, coeffs)
temp2 <- do.call(rbind,pvals)
temp[temp2>psig] = NA
mod_df2 <- data.frame(temp, row.names = colnames(Expr.dat.log))
mod_df2$max_region_coeff = pmax(abs(mod_df2$regionCrossAreal_MTG), 
                                abs(mod_df2$regionCrossAreal_V1), na.rm = TRUE) 
mod_df2 <- mod_df2[order(-mod_df2$max_region_coeff),]
library(tidyverse)
mod_df3 <- mod_df2 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
ggplot(mod_df3, aes(x = colname, y = rowname, fill = value)) +
  geom_tile() +
  labs(title = "Heatmap",
       x = "Coeff",
       y = "Gene")
# Try clustering heatmap
# Sort by descending for max V1 and MTG coeff, take top 20
# Sort by descending for mean abs 

sum(unlist(regional))

# Sanity check expression
foo <- merge(Expr.dat.df, annoAll, by.x = 'row.names', by.y = 'sample_id')
log2fc = log2(mean(foo['CACNB2'][foo['ann_source']=='CrossAreal_M1']) + 1) -
                log2(mean(foo['CACNB2'][foo['ann_source']=='CrossAreal_MTG']) + 1)
# Note: this only checks if normalization by overall expression level is turned off

library(patchwork)
subclass_short = c("D1-Matrix", "D2-Matrix", "D1-Striosome", "D2-Striosome", "D1-ShellOT",
             "D2-ShellOT", "D1D2-Hybrid", "D2-Hybrid-MCHR2", "D1-NUDAP", 
             "PVALB-COL19A1-ST18", "SST_Chodl", "CCK-FBXL7", "CHAT", "CCK-VIP-TAC3", 
              "LHX6-TAC3-PLPP4", "TAC3-LHX8-PLPP4") 
genes = c("SCN1A", "SCN9A", "SCN4B", "KCNQ5", "KCNAB1", "KCNB2", "KCNC2", 
             "KCND2", "KCNH1", "KCNH7", "KCNJ3", "KCNK2", "KCNT2", "KCNIP4", "HCN1", 
             "CACNA1A", "CACNA1D", "ANO10", "ANO3", "ANO4", "ANO6")
genes = c("SCN9A", "SCN4B", "KCNQ5", "KCNAB1", "KCNC2", 
             "KCND2", "KCNH1", "KCNJ3", "KCNK2", "KCNN3", "KCNT2", "KCNIP4",
             "CACNA2D3", "HCN1") 

level = 'level3.subclass_label'
anno_type = anno[[level]]
anno_type_color = anno[[level %>% str_replace("label", "color")]]
keepinds = is.element(anno_type, subclass_short)
Expr.dat <- Expr.dat[keepinds,]  # May be other dimension if running from contiguous
#Expr.dat <- t(Expr.dat)
anno_type <- anno_type[keepinds]
anno_type_color <-anno_type_color[keepinds]
#genes = c("KCNC2", "KCNQ5", "KCNIP4", "CACNA2D3")
colz = anno_type_color[match(subclass_short, anno_type)]

#anno_type <- str_wrap(anno_type, width = 5)
p = list()
for (i in 1:length(genes)){
    gene = genes[i]
    #scale_fill_brewer()
    df <- data.frame(expr = Expr.dat[,gene], subclass = anno_type, subclass_color = anno_type_color)
    df <- mutate(df, expr_log = log2(expr + 1)) 
    df$subclass <- factor(df$subclass, levels = rev(subclass_short))
    # Basic violin plot
    #p[[i]] <- ggplot(df, aes(x=subclass, y=expr_log)) + geom_violin(bw = 0.3) + geom_boxplot(width=0.1) 
    p[[i]] <- ggplot(df, aes(x=subclass, y=expr_log, fill = subclass)) + geom_violin(bw = 0.3, color = 'darkgray', show.legend= FALSE) + scale_fill_manual(values = rev(colz))
    p[[i]] <- p[[i]] + ylab(gene) + xlab(NULL) +
    theme(axis.title.y = element_text(angle=0, size = 36), axis.title.x = element_text(size = 36), axis.text.y = element_text(angle = 0, size = 36)) +
     scale_y_discrete(position = "right") + coord_flip()
    if (i > 1) {
      p[[i]] <- p[[i]] + theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
    }
    
}   
png(file.path(mappingFolder, "/DEG/violins.png"), width = 3200, height = 1100)
p[[1]] | p[[2]] | p[[3]] | p[[4]] | p[[5]] | p[[6]] | p[[7]] | p[[8]] | p[[9]] | p[[10]] | p[[11]] | p[[12]] | p[[13]] | p[[14]]
dev.off()
#+ scale_x_discrete(limits = subclass_short) + scale_y_discrete(position = "top") 


#grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], nrow=4)
#do.call(grid.arrange,c(p, ncol=4))

# PRETTY UP COLORS
# GET DOCALL TO WORK
# COMPARE MSN vs non-MSNlibrary

#BiocManager::install("monocle")
library(monocle)

#Expr.dat <- AIT.anndata$layers['counts']

Expr.sub <- Expr.dat[,Markers_short]
Expr.sub <- Expr.sub[is.element(anno$level3.subclass_label,subclass_short),]
Expr.sp <- t(as.matrix(Expr.sub))
write.csv(Expr.sp, file.path(mappingFolder, 'DEG/Expr_short.csv'))
Expr.sp <- as(Expr.sp, "sparseMatrix") 

fd <- data.frame(gene_short_name=rownames(Expr.sp))
anno <- as.data.frame(anno[is.element(anno$level3.subclass_label,subclass_short),])
write.csv(Expr.sp, file.path(mappingFolder, 'DEG/Anno_short.csv'))
rownames(anno) = colnames(Expr.sp) 
rownames(fd) = rownames(Expr.sp)
cds <- newCellDataSet(Expr.sp,
                         phenoData = new("AnnotatedDataFrame", data = anno),
                         featureData = new("AnnotatedDataFrame", data = fd))


#cds_subset <- cds[row.names(subset(Expr.sp,
#                 gene_short_name %in% allMarkers)),]
plot_genes_violin(cds, min_expr=0.1)

# group_cells_by=anno["level3.subclass_Tree"], 

#Need Monocle3 - hard to install

#########

#Expr = read.csv("/Users/xiaoping.liu/celltypes/NHP_BG_anal/DEG/Expr_short.csv")
#brain.data     <- Expr  # This is already subsetted for cells and
#genes, should be features are rows, cells are columns
#Anno = read.csv("/Users/xiaoping.liu/celltypes/NHP_BG_anal/DEG/Anno_short.csv")
 # rows are cells, columns are metadata fields
#brain.metadata <- data.frame(subclass = Anno$level3.subclass_label)
#rownames(brain.metadata) <- colnames(brain.data)


detach("package:scrattch.bigcat", unload = TRUE)   # detach all scrattch packages on right panel
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("multtest")

source("https://z.umn.edu/archived-seurat")    # Get Seurat V2


## Construct data set lists
brain.data <- Expr.dat
brain.metadata <- anno

brain <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)
# Should you transform counts?
Idents(brain) <- brain.metadata$level3.subclass_label
VlnPlot(object = brain, features = 'counts', idents = subclass_short)