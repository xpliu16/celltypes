## Load scrattch.mapping
library(scrattch.mapping)

## Load in example count data
library(tasic2016data)

## Compute log CPM
query.data = tasic_2016_counts
query.data = logCPM(query.data)

print("Complete")
