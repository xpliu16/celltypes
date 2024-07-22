## Load scrattch.mapping
library(scrattch.mapping)

## Load in example count data
library(tasic2016data)

## Compute log CPM
query.data = tasic_2016_counts
query.data = logCPM(query.data)
annotations <- tasic_2016_anno
# Variable renaming
clusters            <- annotations$primary_type_label[match(1:49,annotations$primary_type_id)]
annotations$cluster <- factor(annotations$primary_type_label, levels=clusters)
annotations$class   <- annotations$broad_type
annotations$class[!is.element(annotations$class,c("GABA-ergic Neuron","Glutamatergic Neuron"))] = "Non-neuronal"


## Specify which taxonomies to map against.
taxonomies = c("/home/docker/tasic_2016")

## Map data against all taxonomies
mapping.anno = list()

mappingFolder = "/home/docker/mapping/"
# mappingFolder = "/home/xiaoping.liu/scrattch/mapping1"

for(taxonomy in taxonomies){
    ## Load the shiny taxonomy into a standard object for mapping.
    AIT.anndata = loadTaxonomy(refFolder = taxonomy)
    ## Map! Returns a data.frame with mapping results.
    mapping.anno[[taxonomy]] = taxonomy_mapping(AIT.anndata=AIT.anndata, 
                                                query.data=query.data, 
                                                label.cols="cluster_label", ## Which obs in AIT.anndata contain annotations to map. E.g. "class", "subclass", etc.
                                                corr.map=TRUE, 
                                                tree.map=TRUE, 
                                                seurat.map=TRUE)
    buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                          mappingFolder  = mappingFolder,
                          query.data     = query.data,  
                          query.metadata = annotations,
                          query.mapping  = mapping.anno[[taxonomy]],
                          doPatchseqQC   = FALSE  # Set to FALSE if not needed or if writePatchseqQCmarkers was not applied in reference generation
    )    
}
