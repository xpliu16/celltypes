Mapmycells

conda activate mapmycells

# Specific to non-neuronal
python -m cell_type_mapper.cli.from_specified_markers \
--query_path /home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116/20240621_RSC-204-366_query.h5ad \
--type_assignment.normalization log2CPM \
--precomputed_stats.path /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/HANN_mapping/Taxonomies/NHP_BG_AIT11.6_neurons/BG_AIT116_neurons_precompute_stats.h5 \
--query_markers.serialized_lookup /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/HANN_mapping/Taxonomies/NHP_BG_AIT11.6_neurons/BG_AIT116_neurons_query_markers.json \
--extended_result_path /home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116/20240621_RSC-204-366_neurons_results.json \
--type_assignment.n_processors 16 \
--type_assignment.chunk_size 1000








****

python -m cell_type_mapper.cli.from_specified_markers \
--query_path /home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116/20240621_RSC-204-366_query.h5ad \
--type_assignment.normalization log2CPM \
--precomputed_stats.path /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/HANN_mapping/Taxonomies/NHP_BG_AIT11.6/BG_AIT116_precompute_stats.h5 \
--query_markers.serialized_lookup /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/HANN_mapping/Taxonomies/NHP_BG_AIT11.6/BG_AIT116_query_markers.json \
--extended_result_path /home/xiaoping.liu/scrattch/mapping/NHP_BG_AIT_116/20240621_RSC-204-366_results.json \
--type_assignment.n_processors 16 \
--type_assignment.chunk_size 1000 

