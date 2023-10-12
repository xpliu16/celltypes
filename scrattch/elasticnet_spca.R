dataFolder = "/Users/xiaoping.liu/celltypes/NHP_BG_anal"

#install.packages('elasticnet')
#install.packages('arrow')
library(feather)
library(elasticnet)

df <- arrow::read_feather(file.path(dataFolder, "X_std_BG.feather"))
#df <- arrow::read_feather(file.path(dataFolder, "X_std.feather"))
ncomp = 92
#para = 600 for current results
sparse.pca.result  <-  spca(df, K = ncomp, para=rep(300, ncomp), type = "predictor", sparse=c("penalty"), lambda=1e-6)
# try setting para to about 5 for each and using varnum - then plot as you change ncomp total var explained

dim(sparse.pca.result$loadings)
length(which(sparse.pca.result$loadings[,1]!=0,arr.ind = T))  # 1->87,  10-> 47,  100-> 13, 1000 -> 9
sparse.pca.result$pev

lst <- sort(sparse.pca.result$pev, index.return=TRUE, decreasing=TRUE)
lst$ix[1:20]

i = 2
x = sparse.pca.result$loadings[,i]    # Scan through components
x[x!=0]

arrow::write_feather(
  as.data.frame(sparse.pca.result$loadings[,1:29]),
  file.path(dataFolder, 'spca_comps.feather')
)

