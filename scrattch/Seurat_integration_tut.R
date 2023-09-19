install.packages("Seurat")
library(Seurat)
#install.packages("SeuratData")
#install.packages("devtools")
#library(devtools)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(patchwork)

getOption('timeout')
options(timeout=200)
InstallData("ifnb")

# Tutorial per: https://satijalab.org/seurat/articles/integration_introduction.html
# load dataset
LoadData("ifnb")
#ifnb$stim

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")
# length(ifnb.list)
