## Load libraries
library(scrattch.taxonomy)
library(scrattch.mapping)
library(scrattch.patchseq)
library(reticulate) # For hierarchical mapping
cell_type_mapper <- import("cell_type_mapper") # For hierarchical mapping

## Specify which reference taxonomy to map against.
## -- Replace folder and file name with correct location
taxonomyDir = getwd() 