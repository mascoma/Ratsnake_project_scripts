#########comparing entire dataset trees
library(ape)
library(outliers)
bsfile1 <- list.files(path="/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/BSbs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
ppfile1 <- list.files(path="/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/pp1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
readtable2<-function (file) {read.table(file,header=T)}
bsgenes<-lapply(bsfile1, function(x) readtable2(x))
ppgenes<-lapply(ppfile1, function(x) readtable2(x))


 
plotdist<-function(files, a, b, name, color){
  matrix<-files[[a]]
  vector<-matrix[,b] 
  if (a==1 && name == "bsgenes"){
    plot(density(vector), xlim=c(0, 1.5), ylim=c(0, 5), xlab="Normalized among-loci tree distance" , col=color, main = NA)
  }
  else {
    lines(density(vector), col=color)
  }
  
}
 
tiff(file="/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/allgenedist_tt_yule.tiff",width=500, height=500)
 
### ttiwith yule model

for (a in 1:300){
  plotdist(bsgenes,a,11,"bsgenes","dodgerblue4")
  plotdist(ppgenes,a, 11, "ppgenes","coral3")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
#lines(density(mlgenes[,14]), lwd = 3, col="chartreuse4")
#lines(density(ppcongenes[,5]),lwd = 3, col="gray25")
legend("topright", legend = c("BS", "PP"),lty=c(1,1,1),lwd= 2, col=c("dodgerblue4", "coral3") )
dev.off()

t.test(bsgenes[[1]][,11], ppgenes[[1]][,11])
###
##Welch Two Sample t-test

###data:  bsgenes[[1]][, 11] and ppgenes[[1]][, 11]
##t = 76.137, df = 91565, p-value < 2.2e-16
##alternative hypothesis: true difference in means is not equal to 0
##95 percent confidence interval:
##  0.08598086 0.09052462
##sample estimates:
##  mean of x mean of y 
##0.5472814 0.4590287 





