mltree <-read.table("E:/documents-export-2015-09-15/forgraph/mldist.txt",header=T)
BItree <-read.table("E:/documents-export-2015-09-15/forgraph/BIdist.txt", header = T)  

mltree_noclk <-read.table("E:/documents-export-2015-09-15/mltree_dists_subset.txt",header=T)
BItree_noclk <-read.table("E:/documents-export-2015-09-15/BItree_dists_subset.txt", header = T)  

mltree_clk <-read.table("E:/documents-export-2015-09-15/ml_dists_rest.txt",header=T)
BItree_clk <-read.table("E:/documents-export-2015-09-15/BI_dists_rest.txt", header = T)  


head(mltree)
mltt<-mltree[[14]]
quantile(mltt)
mean(mltt)
median(mltt)
sd(mltt)
range(mltt)
plot(density(mltt))
 


contt<-BItree[[5]]
quantile(contt)
mean(contt)
sd(contt)
median(contt)
range(contt)
plot(density(contt))
 

mlnoclk<-mltree_noclk[[5]]
quantile(mlnoclk)
mean(mlnoclk)
sd(mlnoclk)
median(mlnoclk)
range(mlnoclk)
plot(density(mlnoclk))



BInoclk<-BItree_noclk[[5]]
quantile(BInoclk)
mean(BInoclk)
sd(BInoclk)
median(BInoclk)
range(BInoclk)
plot(density(BInoclk))

mlclk<-mltree_clk[[5]]
quantile(mlclk)
mean(mlclk)
sd(mlclk)
median(mlclk)
range(mlclk)
plot(density(mlclk))



BIclk<-BItree_clk[[5]]
quantile(BIclk)
mean(BIclk)
sd(BIclk)
median(BIclk)
range(BIclk)
plot(density(BIclk))



plot(density(BInoclk), xlim=c(0, 1.5), ylim=c(0, 5), xlab="normalized tree distance" , col="blue", main = NA)
lines(density(BIclk), lwd = 3, col="chartreuse4")
lines(density(contt), lwd = 3, col="chartreuse4")
