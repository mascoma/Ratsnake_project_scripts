n = 100
g = 6 
set.seed(g)
d <- data.frame(x = unlist(lapply(1:g, function(i) rnorm(n/g, runif(1)*i^2))), 
                y = unlist(lapply(1:g, function(i) rnorm(n/g, runif(1)*i^2))))

require(vegan)
mydata0<-read.table("/Users/Xin/Documents/Ratsnake_project/output/mcmctree/genelist.txt",header=T,row.names=1,sep='\t')
mydata <- na.omit(mydata0)
d = mydata
plot(d)
fit <- cascadeKM(scale(d, center = TRUE,  scale = TRUE), 1, 200, iter = 1000)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

library(mclust)
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20
d_clust <- Mclust(as.matrix(d), G=1:200)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 4 clusters
plot(d_clust)


library(cluster)
clusGap(d, kmeans, 60, B = 100, verbose = interactive())


library(NbClust)
nb <- NbClust(d, diss="NULL", distance = "euclidean", 
              min.nc=2, max.nc=15, method = "kmeans", 
              index = "alllong", alphaBeale = 0.1)
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))