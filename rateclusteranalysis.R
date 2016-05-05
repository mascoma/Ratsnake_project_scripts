mydata0<-read.table("/Users/Xin/Documents/Ratsnake_project/output/mcmctree/group_method4/genelist_reduced.csv",header=T,row.names=1,sep=',')
mydata <- na.omit(mydata0) # listwise deletion of missing
#mydata <- scale(mydata) # standardize variables

## method 1 K-cluster
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:60) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:60, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
fit15=kmeans(mydata, 15)
fit20=kmeans(mydata, 20)
fit10=kmeans(mydata, 10)
aggregate(mydata,by=list(fit20$cluster),FUN=mean)
mydata20 <- data.frame(mydata, fit20$cluster)
write.table(mydata20, file = "/Users/Xin/Documents/Ratsnake_project/output/mcmctree/group4.txt", sep= "\t")
# K-Means Cluster Analysis
fit <- kmeans(mydata, 15) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

## method 2
# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram
groups <- cutree(fit, k=25) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=30, border="red")
# Ward Hierarchical Clustering with Bootstrapped p values
write.table(groups, file = "group2.txt", sep= "\t")

library(pvclust)
fit <- pvclust(mydata, method.hclust="ward.D2",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# method 3
# Model Based Clustering


# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 15)
fit1 <- kmeans(mydata, 9)
fit2 <- kmeans(mydata, 15)
# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster)
# comparing 2 cluster solutions
library(fpc)
cluster.stats(d, fit1$cluster, fit2$cluster)