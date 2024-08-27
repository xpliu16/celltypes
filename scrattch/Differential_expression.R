#refFolderList <- list("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_AIT_116")
refFolderList <- list("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_MTG",
"/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_V1", 
"/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_M1", 
"/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_DLPFC")
#refFolderList <- list("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/CrossAreal_V1") # For within region cluster comparison
mappingFolder <- "/home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116/"

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
  anndataList[[annName]] <- loadTaxonomy(refFolder, 'AI_taxonomy.h5ad')
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
#subset = list('CrossArea_subclass_label', c('L2/3 IT', 'L4 IT'))
subset = list('CrossArea_subclass_label', c('Sst'))
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
                        #title = 'L23IT:MTG_vs_M1',
                        title = 'Sst'
                        splitby = 'CrossArea_subclass_label')) 

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
#figpath = "/home/xiaoping.liu/Desktop/wDLPFC_M1_L23IT6ref"
#figpath = "/home/xiaoping.liu/Desktop/wDLPFC_MTG_L23IT6ref_nointeraction"
figpath = "/home/xiaoping.liu/Desktop/wDLPFC_SST"
dir.create(figpath)

# Cell composition bar plots
df_tallies = data.frame()
for (ann in names(anndataList)) {
  AIT.anndata = anndataList[[ann]]
  tallies <- AIT.anndata$obs[AIT.anndata$obs$CrossArea_subclass_label %in% c ('L2/3 IT', 'L4 IT'),] %>% group_by (CrossArea_cluster_label) %>% tally()
  tallies$n <- tallies$n / sum(tallies$n)
  area = gsub("CrossAreal_","",ann)
  df_temp <- data.frame(x = area, y = tallies$n, Cluster = as.character(tallies$CrossArea_cluster_label))
  df_tallies <- rbind(df_tallies, df_temp)
}
df_tallies <- df_tallies[order(df_tallies$Cluster), ]
png(file.path(figpath, "L23L4IT_composition.png"), width=10, height=5.9, units='in', res=600)
ggplot(df_tallies, aes(x = x, y = y, fill = Cluster)) + 
  geom_bar(stat="identity") + ylab("Fraction") + xlab('') + 
  theme_gray(base_size = 24) + scale_fill_brewer(palette='Set3')
dev.off()

library(lme4)
library (lmerTest) 
library(dplyr)
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
colnames(model_data) = c(colnames(Expr.dat.log), 'region', 'cluster', 'donor')

# Set reference classes
model_data$region <- relevel(factor(model_data$region), ref = "CrossAreal_MTG")
model_data$cluster <- relevel(model_data$cluster, ref = "L2/3 IT_6")
model_data$cluster <- relevel(model_data$cluster, ref = "Sst_1")

# Subsampling cells
reg_counts <- model_data %>% group_by(region) %>% tally()
target_n = min(reg_counts$n)
tallies <- model_data %>% group_by(region, cluster, donor) %>% tally()
tallies <- tallies[order(tallies$n, decreasing=TRUE),]
#counts_M1 <- tallies[tallies$region=='CrossAreal_M1',]
trimmed_data <- model_data

for (reg in unique(reg_counts$region)){
  n_remove = reg_counts$n[reg_counts$region==reg]-target_n
  #print(reg)
  #print(paste0('n_remove ', n_remove))
  to_remove = n_remove
  counts_reg <- tallies[tallies$region==reg,]
  #print(counts_reg)
  n = 1
  while (to_remove > 0){
    if (n < dim(counts_reg)[1]){
      diff1_2 = counts_reg$n[n]-counts_reg$n[n+1]
      print(paste0('diff1_2 ', diff1_2))
      print(dim(trimmed_data)[1])
      samp_n = min(diff1_2*n, to_remove)}
    else{
      samp_n = to_remove
    }
    print(samp_n)
    #data_top <- trimmed_data %>%    # Not quite right, want to filter on list of tuple values
    #  filter(cluster %in% c(counts_reg$cluster[1:n]) 
    #         & region == reg 
    #         & donor %in% c(counts_reg$donor[1:n]))
    data_top <- trimmed_data [(trimmed_data$region == reg) &
                                !is.na(match(mapply(list,trimmed_data$cluster, trimmed_data$donor, SIMPLIFY=FALSE),
                                mapply(list, counts_reg$cluster[1:n], counts_reg$donor[1:n], SIMPLIFY=FALSE))),] 
    temp <- sample(dim(data_top)[1], samp_n, replace = FALSE, prob = NULL)
    trimmed_data <- trimmed_data [!(rownames(trimmed_data) %in% rownames(data_top[temp,])) ,]
    
    #tallies2 <- trimmed_data %>% group_by(region, cluster, donor) %>% tally()
    #tallies2 <- tallies2[order(tallies2$n, decreasing=TRUE),]
    #counts_trimmed <- tallies2[tallies2$region==reg,]
    #print(counts_trimmed)
    
    n <- n+1
    to_remove <- to_remove - samp_n
    print(to_remove)
  }
}

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

# Remove any groups with fewer than 20 samples
tallies2 <- trimmed_data %>% group_by(region, cluster) %>% tally()
tallies3 <- tallies2[tallies2$n<20,]
trimmed_data <- trimmed_data [is.na(match(mapply(list,trimmed_data$cluster, trimmed_data$region, SIMPLIFY=FALSE),
                             mapply(list, tallies3$cluster, tallies3$region, SIMPLIFY=FALSE))),] 
model_data = trimmed_data

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

#for (gene in c('KCNQ5')){
for (gene in colnames(Expr.dat.log)){
  gene <- as.symbol(gene)
  #mod<-eval(bquote(glm(.(gene) ~ region*cluster, family = poisson, data = model_data)))
  #mod<-eval(bquote(glmer(.(gene) ~ region*cluster + (1|donor), family = poisson, data = model_data, nAGQ = 0))) # These pvalues are Wald tests of the model with and without each individual parameter (but not that whole categorical variable)
  mod<-eval(bquote(lmer(.(gene) ~ region + cluster + (1|donor), data = model_data)))
  #mod<-eval(bquote(lme(.(gene) ~ region*cluster, random=~1|donor, data = model_data)))  # this can't handle rank deficiency from not all combinations of crossed variable having samples
  summary(mod)  
  #mod$coefficients
  #mod$residuals
  coeffs <- c(coeffs, list(coef(summary(mod))[,1]))  # coefficient estimates
  pval_i<-coef(summary(mod))[,5]  # column 4 for glm(), column 5 for lmer()
  pvals <- c(pvals, list(pval_i))  # pvals
  
  regional <- c(regional, 
                any((abs(coef(summary(mod))[c('regionCrossAreal_M1',
                                             'regionCrossAreal_V1', 
                                             'regionCrossAreal_DLPFC'),1]) > coeff_thresh)
                    & (p.adjust(pval_i[c('regionCrossAreal_M1',
                                           'regionCrossAreal_V1',
                                         'regionCrossAreal_DLPFC')], 
                                method = "bonferroni") < psig)))
  #terms = rownames(summary(mod)$coefficients)
  #terms_sub = c(terms[grepl("^cluster", terms)])
  cluster <- c(cluster, 
                any((abs(coef(summary(mod))[c('clusterL2/3 IT_5', 
                                             'clusterL2/3 IT_1',
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
                                            'clusterL2/3 IT_1',
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
  
  
  #terms_sub = c(terms[grepl("regionCrossAreal_M1:cluster", terms)],
  #              terms[grepl("regionCrossAreal_V1:cluster", terms)],
  #              terms[grepl("regionCrossAreal_DLPFC:cluster", terms)])
  #regxclust <- c(regxclust,
  #               any((abs(coef(summary(mod))[c(terms_sub),1]) > coeff_thresh)
  #                   & (p.adjust(pval_i[c(terms_sub)], method = "bonferroni") < psig)))
  
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
  donor_ef <- c(donor_ef, ranef(mod))
}  

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
  pval_i<-coef(summary(mod))[,5]  # column 4 for glm(), column 5 for lmer()
  pvals <- c(pvals, list(pval_i))  # pvals
  
  regional <- c(regional, 
                any((abs(coef(summary(mod))[c('regionCrossAreal_MTG',
                                              'regionCrossAreal_V1', 
                                              'regionCrossAreal_DLPFC'),1]) > coeff_thresh)
                    & (p.adjust(pval_i[c('regionCrossAreal_MTG',
                                         'regionCrossAreal_V1',
                                         'regionCrossAreal_DLPFC')], 
                                method = "bonferroni") < psig)))
  v1_spec = c(3, 6, 14, 20, 34, 15)
  terms1 = paste0('clusterSst_', setdiff(2:37,v1_spec))
  terms2 = paste0('clusterSst_', v1_spec)
  terms = c(terms1, terms2)

  cluster <- c(cluster, 
               any((abs(coef(summary(mod))[c(terms),1]) > coeff_thresh)
                   & (p.adjust(pval_i[c(terms)],
                               method = "bonferroni") < psig)))
  reg_clust <- c(reg_clust, 
                 any((abs(coef(summary(mod))[c(terms),1]) > coeff_thresh)
                     & (p.adjust(pval_i[c(terms)],
                                 method = "bonferroni") < psig)))
  
  terms = rownames(summary(mod)$coefficients)
  terms_sub = c(terms[grepl("regionCrossAreal_MTG:cluster", terms)],
                terms[grepl("regionCrossAreal_V1:cluster", terms)],
                terms[grepl("regionCrossAreal_DLPFC:cluster", terms)])
  regxclust <- c(regxclust,
                 any((abs(coef(summary(mod))[c(terms_sub),1]) > coeff_thresh)
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
  donor_ef <- c(donor_ef, ranef(mod))
}  

mod_df <- data.frame(row.names = colnames(Expr.dat.log), 
                     regional = unlist(regional), 
                     cluster = unlist(cluster), 
                     reg_cluster = unlist(reg_clust), 
                     #regxclust = unlist(regxclust),
                     Pr_regional = unlist(Pr_regional),
                     Pr_cluster = unlist(Pr_cluster))
                     #chisq_regional = unlist(chisq_regional),
                     #chisq_cluster = unlist(chisq_cluster))

donor = do.call(cbind, donor_ef)
donor_df = data.frame(t(donor))
p1 <- ggplot(donor_df, aes(x=X1)) + geom_histogram() + scale_x_continuous(limits = c(-1, 1))
p2 <- ggplot(donor_df, aes(x=X2)) + geom_histogram() + scale_x_continuous(limits = c(-1, 1))
p3 <- ggplot(donor_df, aes(x=X3)) + geom_histogram() + scale_x_continuous(limits = c(-1, 1))
p4 <- ggplot(donor_df, aes(x=X4)) + geom_histogram() + scale_x_continuous(limits = c(-1, 1))
p5 <- ggplot(donor_df, aes(x=X5)) + geom_histogram() + scale_x_continuous(limits = c(-1, 1))
png(file.path(figpath, "donor_histogram.png"), width=5, height=4, units='in', res=150)
plot_grid(p1, p2, p3, p4, p5, labels=c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5"), ncol = 1, nrow = 5)
dev.off()

library(VennDiagram)
set0 <- colnames(Expr.dat.log)
set1 <- colnames(Expr.dat.log)[unlist(regional)]
set2 <- colnames(Expr.dat.log)[unlist(cluster)]
set3 <- colnames(Expr.dat.log)[unlist(reg_clust)]
set4 <- colnames(Expr.dat.log)[unlist(regxclust)]
set5 <- colnames(Expr.dat.log)[unlist(Pr_regional)]
set6 <- colnames(Expr.dat.log)[unlist(Pr_cluster)]
myCol <- brewer.pal(3, "Set3")
v<-venn.diagram(
  x = list(set0, set1, set2),
  category.names = c("All" , "Region_coef" , "Cluster_coef"),
  filename = NULL,
  output=TRUE,
  fill=myCol,
  lty = 'blank',
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer",
)
png(file.path(figpath, "venn_diagram1.png"), width=5, height=4, units='in', res=150)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.9, "npc")))
grid.draw(v)
dev.off()

myCol <- brewer.pal(3, "Set3")
v <-venn.diagram(
  x = list(set1, set3, set4),
  category.names = c("Region" , "Region-specialized cluster", "Region x cluster"),
  filename = NULL,
  output=TRUE,
  fill=myCol,
  lty = 'blank',
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)
png(file.path(figpath, "venn_diagram2.png"), width=5, height=4, units='in', res=150)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.9, "npc")))
grid.draw(v)
dev.off()

myCol <- brewer.pal(4, "Set3")
v<-venn.diagram(
  x = list(set1, set2, set5, set6),
  category.names = c("Region_coef", "Cluster_coef", "Region_ANOVA", "Cluster_ANOVA"),
  filename = NULL,
  output=TRUE,
  fill=myCol,
  lty = 'blank',
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)
png(file.path(figpath, "venn_diagram3.png"), width=5, height=4, units='in', res=150)
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.9, "npc")))
grid.draw(v)
dev.off()

# Count
sum(unlist(Pr_cluster))

library(pheatmap)
mod_df <- cbind(mod_df,t(do.call(cbind, Pr_F)), do.call(rbind,pvals))
# If probabilities are not significant, blank out the coefficient as NA
temp <- do.call(rbind, coeffs)
temp2 <- do.call(rbind,pvals)
temp[temp2>psig] = NA     
temp <-temp[ , -which(colnames(temp) == "(Intercept)")]
#col_order = c(1:13,14,16,18,20,15,17,19)  # before DLPFC
#col_order = c(1:14,30,31,15,18,21,24,16,19,22,25,26,17,20,23,27,28,29)
#temp <- temp[, col_order]

mod_df2 <- data.frame(temp, row.names = colnames(Expr.dat.log))
mod_df2$max_region_coeff = pmax(abs(mod_df2$regionCrossAreal_M1), 
                                abs(mod_df2$regionCrossAreal_V1), 
                                abs(mod_df2$regionCrossAreal_DLPFC), na.rm = TRUE) 
mod_df2 <- mod_df2[order(-mod_df2$max_region_coeff),]

v1_spec = c(3, 6, 14, 20, 34, 15) # SST
terms1 = paste0('clusterSst_', setdiff(2:37,v1_spec))
terms2 = paste0('clusterSst_', v1_spec)
terms = c(terms1, terms2)
#terms = c('clusterL2/3 IT_2',    # Actually, don't reorder for L2/3 IT
#          'clusterL2/3 IT_3',
#          'clusterL2/3 IT_4',
#          'clusterL4 IT_2',
#          'clusterL4 IT_3', 
#          'clusterL4 IT_4',
#          'clusterL4 IT_6',
#          'clusterL2/3 IT_1', 
#          'clusterL2/3 IT_5',
#          'clusterL4 IT_1',
#          'clusterL4 IT_5',
#          'clusterL2/3 IT_3', 
#          'clusterL2/3 IT_4')
mod_df2 = mod_df2[,c(colnames(mod_df2[1:3]),terms, colnames(mod_df2[(4+length(terms)):dim(mod_df2)[2]]))]

# Calculate bias to make 0 equal to white on divergent color scale
b <- log(min(mod_df2[1:30,], na.rm=T)/(min(mod_df2[1:30,], na.rm=T)-max(mod_df2[1:30,], na.rm=T)), base = 0.5)
divergent_palette <- colorRampPalette(c("#1b61e4", "white", "coral"), bias=b)
png(file.path(figpath, "Top_30_regional.png"), width=4.5, height=4.25, units='in', res=150) # was width 5 
pheatmap(mod_df2[1:30,1:(dim(mod_df2)[2]-1)], cluster_rows = FALSE, 
         cluster_cols=FALSE, color = divergent_palette(n=100), fontsize=8,
         angle_col=315)
#grid.newpage()
#grid.draw(h)
dev.off()


mod_df2$max_regclust_coeff = pmax(abs(mod_df2$clusterL2.3.IT_2), 
                                abs(mod_df2$clusterL2.3.IT_3), 
                                abs(mod_df2$clusterL2.3.IT_4),
                                abs(mod_df2$clusterL4.IT_2),
                                abs(mod_df2$clusterL4.IT_3),
                                abs(mod_df2$clusterL4.IT_4),
                                abs(mod_df2$clusterL4.IT_6),na.rm = TRUE) 
mod_df2 <- mod_df2[order(-mod_df2$max_regclust_coeff),]

mod_df2$max_regclust_coeff = pmax(abs(mod_df2$clusterSst_3), 
                                  abs(mod_df2$clusterSst_6), 
                                  abs(mod_df2$clusterSst_14),
                                  abs(mod_df2$clusterSst_15),
                                  abs(mod_df2$clusterSst_20),
                                  abs(mod_df2$clusterSst_34),na.rm = TRUE) 
mod_df2 <- mod_df2[order(-mod_df2$max_regclust_coeff),]


b <- log(min(mod_df2[1:30,], na.rm=T)/(min(mod_df2[1:30,], na.rm=T)-max(mod_df2[1:30,], na.rm=T)), base = 0.5)
divergent_palette <- colorRampPalette(c("blue", "white", "red"), bias=b)
png(file.path(figpath, "Top_30_regclust.png"), width=9, height=4, units='in', res=150)
pheatmap(mod_df2[1:30,1:(dim(mod_df2)[2]-2)], cluster_rows = FALSE, 
         cluster_cols=FALSE, color = divergent_palette(n=100), fontsize=6)
dev.off()

terms = colnames(mod_df2)
terms_sub = c(terms[grepl("regionCrossAreal_MTG.cluster", terms)],
              terms[grepl("regionCrossAreal_V1.cluster", terms)],
              terms[grepl("regionCrossAreal_DLPFC.cluster", terms)])
mod_df2$max_regxclust_coeff <- apply(X=abs(mod_df2[,terms_sub]), MARGIN=1, 
                                     function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
          
#mod_df2$max_regxclust_coeff = pmax(abs(mod_df2$regionCrossAreal_MTG.clusterL2.3.IT_5), 
#                                  abs(mod_df2$regionCrossAreal_MTG.clusterL2.3.IT_1), 
#                                  abs(mod_df2$regionCrossAreal_MTG.clusterL4.IT_1),
#                                  abs(mod_df2$regionCrossAreal_MTG.clusterL4.IT_5),
#                                  #abs(mod_df2$regionCrossAreal_V1.clusterL2.3.IT_5),
#                                  abs(mod_df2$regionCrossAreal_V1.clusterL2.3.IT_1),
#                                  abs(mod_df2$regionCrossAreal_V1.clusterL4.IT_1), 
#                                  #abs(mod_df2$regionCrossAreal_M1.clusterL2.3.IT_5),
#                                  #abs(mod_df2$regionCrossAreal_M1.clusterL2.3.IT_1),
#                                  #abs(mod_df2$regionCrossAreal_M1.clusterL4.IT_1),
#                                  #abs(mod_df2$regionCrossAreal_M1.clusterL4.IT_5),
#                                  #abs(mod_df2$regionCrossAreal_M1.clusterL2.3.IT_3),
#                                  #abs(mod_df2$regionCrossAreal_M1.clusterL2.3.IT_4),
#                                  #abs(mod_df2$regionCrossAreal_DLPFC.clusterL2.3.IT_3),
#                                  abs(mod_df2$regionCrossAreal_DLPFC.clusterL2.3.IT_4),
#                                  abs(mod_df2$regionCrossAreal_DLPFC.clusterL2.3.IT_1),
#                                  abs(mod_df2$regionCrossAreal_DLPFC.clusterL2.3.IT_5),
#                                  abs(mod_df2$regionCrossAreal_DLPFC.clusterL4.IT_1),
#                                  abs(mod_df2$regionCrossAreal_DLPFC.clusterL4.IT_5),
#                                  #abs(mod_df2$regionCrossAreal_DLPFC.clusterL2.3.IT_2),
#                                  #abs(mod_df2$regionCrossAreal_DLPFC.clusterL4.IT_2),
#                                  #abs(mod_df2$regionCrossAreal_DLPFC.clusterL4.IT_3),  
#                                  #abs(mod_df2$regionCrossAreal_DLPFC.clusterL4.IT_6), 
#                                  na.rm = TRUE) 
mod_df2 <- mod_df2[order(-mod_df2$max_regxclust_coeff),]

b <- log(min(mod_df2[1:30,], na.rm=T)/(min(mod_df2[1:30,], na.rm=T)-max(mod_df2[1:30,], na.rm=T)), base = 0.5)
divergent_palette <- colorRampPalette(c("blue", "white", "red"), bias=b)
png(file.path(figpath, "Top_30_regxclust.png"), width=5, height=4, units='in', res=150)
pheatmap(mod_df2[1:30,1:(dim(mod_df2)[2]-3)], cluster_rows = FALSE, 
         cluster_cols=FALSE, color = divergent_palette(n=100), fontsize=6)
dev.off()

library(RColorBrewer)
png(file.path(figpath, "ExprLMER_clust_heatmap.png"), width = 2000, height = 4000)
temp <- do.call(rbind, coeffs)
mod_df2 <- data.frame(temp[,1:dim(temp)[2]], row.names = colnames(Expr.dat.log))

b <- log(min(mod_df2, na.rm=T)/(min(mod_df2, na.rm=T)-max(mod_df2, na.rm=T)), base = 0.5)
divergent_palette <- colorRampPalette(c("blue", "white", "red"), bias=b)
pheatmap(mod_df2, cluster_cols=FALSE, color = divergent_palette(n=100), fontsize=24)
dev.off()

# Sanity check expression
foo <- merge(Expr.dat.df, annoAll, by.x = 'row.names', by.y = 'sample_id')
log2fc = log2(mean(foo['CACNB2'][foo['ann_source']=='CrossAreal_M1']) + 1) -
                log2(mean(foo['CACNB2'][foo['ann_source']=='CrossAreal_MTG']) + 1)
# Note: this only checks if normalization by overall expression level is turned off

library(umap)
mapping.umap <- umap(Expr.dat.log)    # Only passed ion channel genes
#mapping.umap <- umap(Expr.dat.df)
#all(rownames(Expr.dat.df) == rownames(Expr.dat.log))  # TRUE
layout <- mapping.umap$layout
#save(layout, file="/home/xiaoping.liu/Desktop/L23IT_L4IT_UMAP_allgenes.Rdata")
save(mapping.umap, file=file.path(figpath,"L23IT_L4IT_UMAP_gene_subset.Rdata"))
# Shuffle rows so some categories don't end up all on top in scatterplot
randorder = sample(dim(layout)[1], dim(layout)[1], replace = FALSE)
layout <- layout[randorder,]
annoAll_shuff <- annoAll[randorder,]
Expr.dat.shuff_logCPM <- log2(Expr.dat.cpm+1)[randorder,]
colnames(Expr.dat.shuff_logCPM) <- paste0(colnames(Expr.dat.shuff_logCPM),'_logCPM')
Expr.dat.shuff <- Expr.dat.df[randorder,]    # raw counts

#umap_gene_list = c('KCNIP4', 'KCND3', 'KCND2', 'DPP10', 'KCNMB2', 'KCNQ5', 
#                   'CACNG4', 'SCN3B', 'KCNH7')
umap_gene_list = c('KCNH8', 'KCNIP1', 'KCNMB2', 'SCN3B', 'CACNG3', 
                   'ANO3', 'CNGB1', 'KCNH7')
umap_gene_list_log <- paste(umap_gene_list, "logCPM", sep='_')
umap_df <- data.frame(layout,
                      cluster = annoAll_shuff$CrossArea_cluster_label, 
                      subclass = annoAll_shuff$CrossArea_subclass_label,
                      region = annoAll_shuff$ann_source,
                      Expr.dat.shuff[,umap_gene_list], 
                      Expr.dat.shuff_logCPM[,umap_gene_list_log])

#umap_df$cluster <- factor(umap_df$cluster, levels=c("L2/3 IT_1", "L2/3 IT_2",   
#                                                    "L2/3 IT_3", "L2/3 IT_4",
#                                                    "L2/3 IT_5", "L2/3 IT_6",  
#                                                    "L4 IT_1", "L4 IT_2",
#                                                    "L4 IT_3", "L4 IT_4",
#                                                    "L4 IT_5", "L4 IT_6"))

png(file.path(figpath,"L23IT_L4IT_UMAP_cluster.png"), width=4.8, height=4, units='in', res=150)
theme_set(theme_grey(base_size = 18)) 
#png(file.path(figpath,"SST_UMAP_cluster2.png"), width=5.4, height=4, units='in', res=150)
p1<-ggplot(umap_df, aes(x=X1, y=X2, color=cluster)) + geom_point(size=0.5) + scale_color_brewer(palette="Set3")    
# TECHNICALLY NOT ENOUGH COLORS FOR SST

#p1 + ggtitle("L2/3 IT and L4 IT cells\n(M1, DLPFC, MTG, and V1) Ion channel genes") +
p1 + ggtitle("SST cells\n(M1, DLPFC, MTG, and V1) Ion channel genes") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "#F2F2F2"), 
              axis.line = element_line(color = "black", linewidth = 0.7),
              axis.text.x=element_blank(), 
              axis.ticks.x=element_blank(), 
              axis.text.y=element_blank(), 
              axis.ticks.y=element_blank())
dev.off()


#png(file.path(figpath,"SST_UMAP_subclass.png"), width=5.4, height=4, units='in', res=150)
png(file.path(figpath,"L23IT_L4IT_UMAP_subclass.png"), width=4.5, height=4, units='in', res=150)
p1<-ggplot(umap_df, aes(x=X1, y=X2, color=subclass)) + geom_point(size=0.5) + scale_color_brewer(palette="Set3")

p1 + ggtitle("L2/3 IT and L4 IT cells\n(M1, DLPFC, MTG, and V1) Ion channel genes") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#F2F2F2"), 
        axis.line = element_line(color = "black", linewidth = 0.7),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
dev.off()

#png(file.path(figpath,"SST_UMAP_region.png"), width=5.4, height=4, units='in', res=150)
png(file.path(figpath,"L23IT_L4IT_UMAP_region.png"), width=5.7, height=4, units='in', res=150)
p1<-ggplot(umap_df, aes(x=X1, y=X2, color=region)) + geom_point(size=0.5) + scale_color_brewer(palette="Set3")

p1 + ggtitle("L2/3 IT and L4 IT cells\n(M1, DLPFC, MTG, and V1) Ion channel genes") +
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#F2F2F2"), 
        axis.line = element_line(color = "black", linewidth = 0.7),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
dev.off()

png(file.path(figpath,"L23IT_L4IT_UMAP_KCNMB2.png"), width=4.8, height=4, units='in', res=150)
p1<-ggplot(umap_df, aes(x=X1, y=X2, color=KCNMB2)) + geom_point(size=0.5)

p1 + ggtitle("L2/3 IT and L4 IT cells\n(M1, MTG, and V1) Ion channel genes") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#F2F2F2"), 
        axis.line = element_line(color = "black"))
dev.off()

p1 <- ggplot(umap_df, aes(x=KCNMB2)) + geom_histogram()
p2 <- ggplot(umap_df, aes(x=KCNMB2_logCPM)) + geom_histogram()
plot_grid(p1, p2, ncol = 2, nrow = 1)

plot_umap_expr <- function(umap_df, gene, path = '/home/xiaoping.liu/Desktop') {
  fn = sprintf("L23IT_L4IT_UMAP_%s.png", gene)
  print(file.path(path, fn))
  png(file.path(path, fn), width=3.5, height=4, units='in', res=150)
  gene <- as.symbol(gene)
  p1<-eval(bquote(ggplot(umap_df, aes(x=X1, y=X2, color=.(gene))) + geom_point(size=0.5)))
  
  print(p1 + ggtitle("L2/3 IT and L4 IT cells\n(M1, DLPFC, MTG, and V1) Ion channel genes") +
    xlab("UMAP1") + ylab("UMAP2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#F2F2F2"), 
          axis.line = element_line(color = "black", linewidth = 0.7),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          legend.position = c(0.82, 0.77),
          legend.text = element_text(size = 14),
          legend.title=element_blank(),
          legend.key.size = unit(0.17, "in")))
  
  dev.off()
}

for (gene in umap_gene_list){ 
  plot_umap_expr (umap_df, gene, figpath)
  plot_umap_expr (umap_df, paste0(gene,'_logCPM'), figpath)
}


umap_df$region <- gsub("^CrossAreal_","",umap_df$region)
png("/home/xiaoping.liu/Desktop/L23IT_L4IT_UMAP_KCNMB2_violin.png", width=3, height=3, units='in', res=150)
p <- ggplot(umap_df, aes(x=region, y=KCNMB2)) + geom_violin(bw = 0.3, color = 'darkblue', show.legend= FALSE) + ylab('KCNMB2') + xlab(NULL)
#  theme(axis.title.y = element_text(angle=0), axis.title.x = element_text(), axis.text.y = element_text(angle = 0, size = 36)) +
#  scale_y_discrete(position = "right") + coord_flip()
p + stat_summary(fun.y=median, geom="point", size=2, color="red")
dev.off()


df_PCA = read.csv("/home/xiaoping.liu/Desktop/SST_PCA.csv")  
vars <- load("/home/xiaoping.liu/Desktop/wDLPFC_SST/L23IT_L4IT_UMAP_gene_subset.Rdata")    # Misnamed, actually SST
# Load Expr_patchseq
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/"
data_fn = "20240321_RSC-122-359_human_patchseq_star2.7"
vars2 <- load(paste0(data_dir, paste0(data_fn, "_cpm.Rdata")))
# Genes are rows, samples are columns
#query.data <- logCPM(cpmR)
query.data <- log2(cpmR +1)
query.data <- query.data[colnames(Expr.dat.log),]
query.data <- t(query.data)
patchseq.umap <- predict(mapping.umap, query.data)
foo <- merge(x = df_PCA, y = patchseq.umap, by.x = 'exp_component_name', by.y = 0) 

library("viridis")
#png(file.path(figpath,"SST_ephys_PC1.png"), width=4.82, height=4, units='in', res=150)
png(file.path(figpath,"SST_ephys_PC1_legnotitle_logplusone.png"), width=6.025, height=5, units='in', res=150)
p1<-ggplot(foo, aes(x=V1, y=V2, color=pca_comp_1)) + geom_point(size=1.2, alpha=0.9) 
#+ scale_color_brewer(palette="viridis")

print(p1 + ggtitle("") +
        xlab("UMAP1") + ylab("UMAP2") +
        #scale_color_gradient(limits = c(-8, 8)) +
        scale_color_viridis_c(limits = c(-8, 8), oob = scales::squish) +
        labs(color='') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "#F2F2F2"), 
              axis.line = element_line(color = "black"),
              legend.spacing.x = unit(8, 'mm')))
dev.off()

png(file.path(figpath,"SST_ephys_PC2.png"), width=5.2, height=4, units='in', res=150)
p1<-ggplot(foo, aes(x=V1, y=V2, color=pca_comp_2)) + geom_point(size=0.5) 
#+ scale_color_brewer(palette="Set3")

print(p1 + ggtitle("") +
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_gradient(limits = c(-8, 8)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "#F2F2F2"), 
              axis.line = element_line(color = "black")))
  
dev.off()


png(file.path(figpath,"SST_ephys_upstroke_adapt_ratio.png"), width=5.2, height=4, units='in', res=150)
p3<-ggplot(foo, aes(x=V1, y=V2, color=upstroke_adapt_ratio)) + geom_point(size=0.5, alpha = 0.9) 
#+ scale_color_brewer(palette="Set3")

print(p3 + ggtitle("") +
        xlab("UMAP1") + ylab("UMAP2") +
        #scale_color_viridis() +
        scale_color_viridis_c(limits = c(0.3579092, 1.05), oob = scales::squish) +
        labs(color='Upstroke Adapt') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "#F2F2F2"), 
              axis.line = element_line(color = "black")))
dev.off()

png(file.path(figpath,"SST_ephys_downstroke_hero.png"), width=5.2, height=4, units='in', res=150)
p1<-ggplot(foo, aes(x=V1, y=V2, color=downstroke_hero)) + geom_point(size=0.5) 

print(p1 + ggtitle("") +
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_viridis_c(limits = c(-330, -50), oob = scales::squish) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "#F2F2F2"), 
              axis.line = element_line(color = "black")))
dev.off()

foo$width_suprathresh_hero_ms = foo$width_suprathresh_hero*1000
png(file.path(figpath,"SST_ephys_width_suprathresh_hero.png"), width=4.84, height=4, units='in', res=150)
p2<-ggplot(foo, aes(x=V1, y=V2, color=width_suprathresh_hero_ms)) + geom_point(size=0.5, alpha=0.9) 

print(p2 + ggtitle("") +
        xlab("UMAP1") + ylab("UMAP2") +
        #scale_color_viridis() +
        scale_color_viridis_c(limits = c(0.26, 1.2), oob = scales::squish) +
        labs(color='Width (ms)') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "#F2F2F2"), 
              axis.line = element_line(color = "black")))
dev.off()

# Does not work
png(file.path(figpath,"SST_ephys_composite.png"), width=14.5, height=4, units='in', res=150)
grid.arrange(p1, p2, p3, ncol = 3, widths = c(1, 1, 1))
dev.off()

library(egg)
png(file.path(figpath,"SST_ephys_composite.png"), width=14.5, height=4, units='in', res=150)
ggarrange(p1,p2,p3, ncol = 3)
dev.off()

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

#PVALB left and right wings
vars<-load(file.path(mappingFolder,'Int_UMAP/PVALB_left_right.Rdata'))
data_fn = "20240520_RSC-204-363_macaque_patchseq_star2.7"
load(file.path(data_dir, paste0(data_fn, "_cpm.Rdata")))
load(file.path(data_dir, paste0(data_fn, "_samp.dat.Rdata")))

Expr.dat <- AIT.anndata$layers['counts']
#inds = AIT.anndata$obs['Subclass_label']=='PVALB-COL19A1-ST18'
#Expr.dat <- Expr.dat[inds,]
Expr.dat<-t(Expr.dat)
Expr.dat_norm <- Expr.dat

ident_vec <- rep("Other", dim(AIT.anndata$obs)[1])
ident_vec[is.element(AIT.anndata$obs$sample_id, PV_l_tax_names)] = "PV_l"
ident_vec[is.element(AIT.anndata$obs$sample_id, PV_r_tax_names)] = "PV_r"

dataBG_all<-Expr.dat_norm
anno_all<-AIT.anndata$obs

#brain.data     <- cbind(dataBG_all[keepGenes,],dataBG_all_PS[keepGenes,])  # Include only genes subsetted above
brain.data <- dataBG_all
brain.metadata <- data.frame(type = anno_all$Subclass_label)
rownames(brain.metadata) <- colnames(brain.data)

## Construct data set lists
brain      <- CreateSeuratObject(counts = brain.data, meta.data = brain.metadata)
Idents(brain) <- ident_vec

de.markers <- FindMarkers(brain, ident.1 = "PV_l", ident.2 = "PV_r")    # If you want to downsample: max.cells.per.ident = max_n_cells
  
save(de.markers, file = file.path(mappingFolder,'PVALB_lr_DE.Rdata'))

FCcutoff = 0.5  
title = 'PVALB-COL19A1-ST18 Left vs. Right wing'

png(file.path(mappingFolder, paste0('PVALB_lr_DEG.png')), width = 1600, height = 1200)
p1 <- EnhancedVolcano(de.markers,
                      lab = rownames(de.markers),
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

# Prevalence
anno<-AIT.anndata$obs
anno<- anno[anno$Subclass_label=='PVALB-COL19A1-ST18',]
#tallies_m <- anno[anno$Sex_label =='M',] %>% group_by (Sex_label, cluster_label) %>% tally()
tallies <- anno %>% group_by (Sex_label, cluster_label) %>% tally()
tallies = data.frame(tallies)
tallies[tallies['Sex_label']=='M','n'] <- tallies[tallies['Sex_label']=='M','n'] / sum(tallies[tallies['Sex_label']=='M','n'])
tallies[tallies['Sex_label']=='F','n'] <- tallies[tallies['Sex_label']=='F','n'] / sum(tallies[tallies['Sex_label']=='F','n'])
tallies = tallies[tallies$Sex_label != "unknown",]

png(file.path(mappingFolder, "PVALB_composition.png"), width=10, height=5.9, units='in', res=600)
ggplot(tallies, aes(x = Sex_label, y = n, fill = cluster_label)) + 
  geom_bar(stat="identity") + ylab("Fraction") + xlab('Sex') + 
  theme_gray(base_size = 24) + scale_fill_brewer(palette='Set3')
dev.off()

