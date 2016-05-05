#########comparing entire dataset trees

bsfile1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/BSbs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
ppfile1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/pp1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
mlgenes <-read.table("mldist.txt",header=T)
ppcongenes <-read.table("condist.txt", header=T)  

#G6file1 <- list.files(path="C:/Users/Xin/Desktop/sumresults/G6bs1000outforgraph", pattern="*.*", full.names=T, recursive=F)

readtable2<-function (file) {read.table(file,header=T)}

bsgenes<-lapply(bsfile1, function(x) readtable2(x))
ppgenes<-lapply(ppfile1, function(x) readtable2(x))
#G3genes<-lapply(G3file1, function(x) readtable2(x))
#G5genes<-lapply(G5file1, function(x) readtable2(x))
#G6genes<-lapply(G6file1, function(x) readtable2(x))

#wmatrix<- G1genes[[1]]
#vector<-matrix[,5]
# plot(density(vector), xlim=c(0, 1.5), ylim=c(0, 3), xlab="normalized tree distance", col="blue")
#file<-G1genes
#names(file)

#######################################
spcmp<-read.csv("specmp.csv")
spmc<-data.frame(spcmp[,6],spcmp[,3])
plot(density(sqrt(spcmp[,3])))
 
t.test(spmc[,2]~spmc[,1], data=spmc)


sprc<-data.frame(spcmp[,6],spcmp[,4])
plot(density((spcmp[,4])))

t.test(sprc[,2]~sprc[,1], data=sprc)

sptt<-data.frame(spcmp[,6],spcmp[,5])
plot(density(sqrt(spcmp[,5])))

t.test(sptt[,2]~sptt[,1], data=sptt)


####anova test for four gene sets

lmlmc<-log(mlmc<-mlgenes[[5]])
plot(density(lmlmc))
qqnorm(lmlmc)


lconmc<-log(ppcongenes[[5]])
plot(density(lconmc))
qqnorm(lconmc)

lbsmc<-log(bsgenes[[10]][[5]])
plot(density(lbsmc))
qqnorm(lbsmc)


lppmc<-log(ppgenes[[11]][[5]])
plot(density(lppmc))
qqnorm(lppmc) 

plot(density(lmlmc),xlim=c(-5,2))
lines(density(lconmc))
lines(density(lbsmc))
lines(density(lppmc))
t.test(lbsmc,lppmc)

fac<-rep(c("ml","co","bs","pp"),each=46056)
v1<-c(lmlmc,lconmc,lbsmc,lppmc)
mcyule<-data.frame(v1,fac)
head(mcyule)
 
result<-aov(v1~fac,data=mcyule)
summary(result)
TukeyHSD(result,conf.level=0.95)
pairwise.t.test(v1, fac, p.adjust="bonferroni")
plot(v1~fac, data=mcyule)


mlrc<-mlgenes[[8]]
plot(density(sqrt(mlrc)))
qqnorm(sqrt(mlrc))
smlrc<-sqrt(mlrc)


conrc<-ppcongenes[[8]]
plot(density(sqrt(conrc)))
qqnorm(sqrt(conrc))
sconrc<-sqrt(conrc)


bsrc<-bsgenes[[10]][[8]]
plot(density(sqrt(bsrc)))
qqnorm(sqrt(bsrc))
sbsrc<-sqrt(bsrc)


pprc<-ppgenes[[10]][[8]]
plot(density(sqrt(pprc)))
qqnorm(sqrt(pprc)) 
spprc<-sqrt(pprc)



plot(density(smlrc), ylim=c(0,35))
lines(density(sconrc))
lines(density(sbsrc))
lines(density(spprc))
# t.test(lbsmc,lppmc)

fac<-rep(c("ml","co","bs","pp"),each=46056)
v1<-c(smlrc,sconrc,sbsrc,spprc)
rcyule<-data.frame(v1,fac)
head(rcyule)

result<-aov(v1~fac,data=rcyule)
summary(result)
TukeyHSD(result,conf.level=0.95)
pairwise.t.test(v1, fac, p.adjust="bonferroni")

plot(v1~fac, data=rcyule)


mltt<-mlgenes[[14]]
plot(density(sqrt(mltt)))
qqnorm(sqrt(mltt))
smltt<-sqrt(mltt)


contt<-ppcongenes[[14]]
plot(density(sqrt(contt)))
qqnorm(sqrt(contt))
scontt<-sqrt(contt)


bstt<-bsgenes[[10]][[11]]
plot(density(sqrt(bstt)))
qqnorm(sqrt(bstt))
sbstt<-sqrt(bstt)


pptt<-ppgenes[[10]][[11]]
plot(density(sqrt(pptt)))
qqnorm(sqrt(pptt)) 
spptt<-sqrt(pptt)



plot(density(smltt), ylim=c(0,35))
lines(density(scontt))
lines(density(sbstt))
lines(density(spptt))
# t.test(lbsmc,lppmc)

fac<-rep(c("ml","co","bs","pp"),each=46056)
v1<-c(smltt,scontt,sbstt,spptt)
ttyule<-data.frame(v1,fac)
head(ttyule)

result<-aov(v1~fac,data=ttyule)
summary(result)
TukeyHSD(result,conf.level=0.95)
pairwise.t.test(v1, fac, p.adjust="bonferroni")

plot(v1~fac, data=ttyule)


plotdist<-function(files, a, b, name, color){
  matrix<-files[[a]]
  vector<-matrix[,b] 
  if (a==1 && name == "bsgenes"){
    plot(density(vector), xlim=c(0, 2), ylim=c(0, 8), xlab="normalized tree distance", col=color)
  }
  else {
    lines(density(vector), col=color)
  }
}
pdf("allgenedist_mc_yule).pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(bsgenes,a,5,"bsgenes","blue")
  plotdist(ppgenes,a, 5, "ppgenes","red")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
lines(density(mlgenes[,5]), col="green")
lines(density(ppcongenes[,5]),col="black")
dev.off()

pdf("allgenedist_mc_unif.pdf")
#### gene distribution mc with unifor model matrix 
for (a in 1:300){
  plotdist(bsgenes,a,6,"bsgenes","blue")
  plotdist(ppgenes,a, 6, "ppgenes","red")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
lines(density(mlgenes[,6]), col="green")
lines(density(ppcongenes[,6]),col="black")
 
dev.off()

pdf("allgenedist_rc_yule.pdf")
### rc with yule model 
for (a in 1:300){
  plotdist(bsgenes,a,8,"bsgenes","blue")
  plotdist(ppgenes,a, 8, "ppgenes","red")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
lines(density(mlgenes[,8]), col="green")
lines(density(ppcongenes[,8]),col="black")
dev.off()

pdf("allgenedist_rc_unif.pdf")

### rc with unif model 
for (a in 1:300){
  plotdist(bsgenes,a,9,"bsgenes","blue")
  plotdist(ppgenes,a, 9, "ppgenes","red")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
lines(density(mlgenes[,9]), col="green")
lines(density(ppcongenes[,9]),col="black")
dev.off()

pdf("allgenedist_tt_yule.pdf")
### ttiwith yule model

for (a in 1:300){
  plotdist(bsgenes,a,11,"bsgenes","blue")
  plotdist(ppgenes,a, 11, "ppgenes","red")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
lines(density(mlgenes[,11]), col="green")
lines(density(ppcongenes[,11]),col="black")
dev.off()

pdf("allgenedist_tt_unif.pdf")
#####tt with unifr model
for (a in 1:300){
  plotdist(bsgenes,a,12,"bsgenes","blue")
  plotdist(ppgenes,a, 12, "ppgenes","red")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
lines(density(mlgenes[,12]), col="green")
lines(density(ppcongenes[,12]),col="black")
dev.off()


###########################################################3
####################################################
#################################################
###############################################
################ compare gene trees vs specids trees

bsmp <- list.files(path="G:/work/treedistribution/treecmp/sumresults/BSgtvssptoutput/mptreedist", pattern="*.txt", full.names=T, recursive=F)
bsRF <- list.files(path="G:/work/treedistribution/treecmp/sumresults/BSgtvssptoutput/RFtreedist", pattern="*.txt", full.names=T, recursive=F)
ppmp <- list.files(path="G:/work/treedistribution/treecmp/sumresults/PPgtvssptoutput/mptreedist", pattern="*.txt", full.names=T, recursive=F)
ppRF <- list.files(path="G:/work/treedistribution/treecmp/sumresults/PPgtvssptoutput/RFtreedist", pattern="*.txt", full.names=T, recursive=F)
##G3file2 <- list.files(path="C:/Users/Xin/Desktop/sumresults/G3gtvssptoutput", pattern="*.txt", full.names=T, recursive=F)
#G5file2 <- list.files(path="C:/Users/Xin/Desktop/sumresults/G5gtvssptoutput", pattern="*.txt", full.names=T, recursive=F)
#G6file2 <- list.files(path="C:/Users/Xin/Desktop/sumresults/G6gtvssptoutput", pattern="*.txt", full.names=T, recursive=F)
mlmpt <-read.table("G:/work/treedistribution/treecmp/sumresults/mlgtvssptoutput/mlmpdist.txt",header=T)
mlRFt <-read.table("G:/work/treedistribution/treecmp/sumresults/mlgtvssptoutput/mlgtvsRF.txt",header=T)
conmpt <-read.table("G:/work/treedistribution/treecmp/sumresults/congtvssptoutput/ppgetvsmp.txt", header=T)  
conRFt <-read.table("G:/work/treedistribution/treecmp/sumresults/congtvssptoutput/congtvsRF.txt", header=T)  


readtable2<-function (file) {read.table(file,header=T)}

bsmpt<-lapply(bsmp, function(x) readtable2(x))
bsRFt<-lapply(bsRF, function(x) readtable2(x))
ppmpt<-lapply(ppmp, function(x) readtable2(x))
ppRFt<-lapply(ppRF, function(x) readtable2(x))
#G3gene2<-lapply(G3file2, function(x) readtable2(x))
#G5gene2<-lapply(G5file2, function(x) readtable2(x))
#G6gene2<-lapply(G6file2, function(x) readtable2(x))

#wmatrix<- G1genes[[1]]
#vector<-matrix[,5]
# plot(density(vector), xlim=c(0, 1.5), ylim=c(0, 3), xlab="normalized tree distance", col="blue")
#file<-G1genes
#names(file)
plotdist<-function(files, a, b, name1,name2, color, x,y){
  matrix<-files[[a]]
  vector<-matrix[,b] 
  if (a==1 && name1 == name2){
    plot(density(vector), xlim=c(0, x), ylim=c(0, y), xlab="normalized tree distance", col=color)
  }
  else {
    lines(density(vector), col=color)
  }
}

 

####column 3-mc,6-rc,9-tt
pdf("bsgtvsspt_mc_yule.pdf")
for (a in 1:300){
  plotdist(bsmpt, a, 3, "bsmpt","bsmpt","blue",2,4)
  plotdist(bsRFt, a, 3, "bsmpt", "bsRFt", "red",2,4)
  plotdist(bsgenes, a, 5, "bsmpt","bsgenes","black",2,4)
}

dev.off()

pdf("ppgtvsspt_mc_yule.pdf")
for (a in 1:300){
  plotdist(ppmpt, a, 3, "ppmpt","ppmpt","blue",2,4)
  plotdist(ppRFt, a, 3, "ppmpt", "ppRFt", "red",2,4)
  plotdist(ppgenes, a, 5, "ppmpt","ppgenes","black",2,4)
}

dev.off()


pdf("mlgtvsspt_mc_yule.pdf")

  plot(density(mlmpt[,3]), col="blue", xlim=c(0,2),ylim=c(0,4))
  lines(density(mlRFt[,3]), col="red",xlim=2,ylim=4)
  lines(density(mlgenes[,5]), col="black",xlim=2,ylim=4)
 
dev.off()


pdf("congtvsspt_mc_yule.pdf")
 
  plot(density(conmpt[,3]), col="blue",xlim=c(0,2),ylim=c(0,5))
  lines(density(conRFt[,3]), col="red")
  lines(density(ppcongenes[,5]),col="black")

dev.off()  

###############rc


pdf("bsgtvsspt_rc_yule.pdf")
for (a in 1:300){
  plotdist(bsmpt, a, 6, "bsmpt","bsmpt","blue",1.5,15)
  plotdist(bsRFt, a, 6, "bsmpt", "bsRFt", "red",1.5,15)
  plotdist(bsgenes, a, 8, "bsmpt","bsgenes","black",1.5,15)
}

dev.off()

pdf("ppgtvsspt_rc_yule.pdf")
for (a in 1:300){
  plotdist(ppmpt, a, 6, "ppmpt","ppmpt","blue",1.5,15)
  plotdist(ppRFt, a, 6, "ppmpt", "ppRFt", "red",1.5,15)
  plotdist(ppgenes, a, 8, "ppmpt","ppgenes","black",1.5,15)
}

dev.off()


pdf("mlgtvsspt_rc_yule.pdf")

plot(density(mlmpt[,6]), col="blue", xlim=c(0,1.5),ylim=c(0,15))
lines(density(mlRFt[,6]), col="red")
lines(density(mlgenes[,8]), col="black")

dev.off()


pdf("congtvsspt_rc_yule.pdf")

plot(density(conmpt[,6]), col="blue",xlim=c(0,1.5),ylim=c(0,15))
lines(density(conRFt[,6]), col="red")
lines(density(ppcongenes[,8]),col="black")

dev.off()  



##################tt

pdf("bsgtvsspt_tt_yule.pdf")
for (a in 1:300){
  plotdist(bsmpt, a, 9, "bsmpt","bsmpt","blue",2,5)
  plotdist(bsRFt, a, 9, "bsmpt", "bsRFt", "red",2,5)
  plotdist(bsgenes, a, 11, "bsmpt","bsgenes","black",2,5)
}

dev.off()

pdf("ppgtvsspt_tt_yule.pdf")
for (a in 1:300){
  plotdist(ppmpt, a, 9, "ppmpt","ppmpt","blue",2,6)
  plotdist(ppRFt, a, 9, "ppmpt", "ppRFt", "red",2,6)
  plotdist(ppgenes, a, 11, "ppmpt","ppgenes","black",2,6)
}

dev.off()


pdf("mlgtvsspt_tt_yule.pdf")

plot(density(mlmpt[,9]), col="blue", xlim=c(0,2),ylim=c(0,6))
lines(density(mlRFt[,9]), col="red")
lines(density(mlgenes[,11]), col="black")

dev.off()


pdf("congtvsspt_tt_yule.pdf")

plot(density(conmpt[,9]), col="blue",xlim=c(0,2),ylim=c(0,6))
lines(density(conRFt[,9]), col="red")
lines(density(ppcongenes[,11]),col="black")

dev.off()  





































#### mc with yule model matrix

head(bsgene2[[1]])
for (a in 1:600){
  plotdist(bsgene2,a,3,"bsgene2","blue")
  plotdist(ppgene2,a, 3, "ppgene2","red")
  #plotdist(G3gene2,a,3,"G3gene2","green")
  #plotdist(G5gene2,a,3,"G5gene2","black")
  #plotdist(G6gene2,a,3,"G6gene2","purple")
}
lines(density(mlgene2[,3]), col="green")
lines(density(ppcongene2[,3]),col="black")
dev.off()

pdf("allgtspt_mc_unif.pdf")
#### mc with unifor model matrix 
for (a in 1:600){
  plotdist(bsgene2,a,4,"bsgene2","blue")
  plotdist(ppgene2,a, 4, "ppgene2","red")
  #plotdist(G3gene2,a,4,"G3gene2","green")
  #plotdist(G5gene2,a,4,"G5gene2","black")
  #plotdist(G6gene2,a,4,"G6gene2","purple")
}
lines(density(mlgene2[,4]), col="green")
lines(density(ppcongene2[,4]),col="black")
dev.off()

pdf("allgtspt_rc_yule.pdf")
### rc with yule model 
for (a in 1:600){
  plotdist(bsgene2,a,6,"bsgene2","blue")
  plotdist(ppgene2,a, 6, "ppgene2","red")
  #plotdist(G3gene2,a,4,"G3gene2","green")
  #plotdist(G5gene2,a,4,"G5gene2","black")
  #plotdist(G6gene2,a,4,"G6gene2","purple")
}
lines(density(mlgene2[,6]), col="green")
lines(density(ppcongene2[,6]),col="black")
dev.off()

pdf("allgtspt_rc_unif.pdf")

### rc with unif model 
for (a in 1:600){
  plotdist(bsgene2,a,7,"bsgene2","blue")
  plotdist(ppgene2,a, 7, "ppgene2","red")
  #plotdist(G3gene2,a,4,"G3gene2","green")
  #plotdist(G5gene2,a,4,"G5gene2","black")
  #plotdist(G6gene2,a,4,"G6gene2","purple")
}
lines(density(mlgene2[,7]), col="green")
lines(density(ppcongene2[,7]),col="black")
dev.off()

pdf("allgtspt_tt_yule.pdf")
### ttiwith yule model

for (a in 1:600){
  plotdist(bsgene2,a,9,"bsgene2","blue")
  plotdist(ppgene2,a, 9, "ppgene2","red")
  #plotdist(G3gene2,a,4,"G3gene2","green")
  #plotdist(G5gene2,a,4,"G5gene2","black")
  #plotdist(G6gene2,a,4,"G6gene2","purple")
}
lines(density(mlgene2[,12]), col="green")
lines(density(ppcongene2[,12]),col="black")
dev.off()

pdf("allgtspt_tt_unif.pdf")
#####tt with unifr model
for (a in 1:600){
  plotdist(bsgene2,a,10,"bsgene2","blue")
  plotdist(ppgene2,a, 10, "ppgene2","red")
  #plotdist(G3gene2,a,4,"G3gene2","green")
  #plotdist(G5gene2,a,4,"G5gene2","black")
  #plotdist(G6gene2,a,4,"G6gene2","purple")
}
lines(density(mlgene2[,13]), col="green")
lines(density(ppcongene2[,13]),col="black")
dev.off()


##########################
#########################


bsfile3 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/bstrees304output", pattern="*.txt", full.names=T, recursive=F)
ppfile3 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/pptrees304output", pattern="*.txt", full.names=T, recursive=F)
readtable2<-function (file) {read.table(file,header=T)}

bsgene3<-lapply(bsfile3, function(x) readtable2(x))
ppgene3<-lapply(ppfile3, function(x) readtable2(x))

plotdist<-function(files, a, b, name, color,x,y){
  matrix<-files[[a]]
  vector<-matrix[,b] 
  if (a==1 && name == "bsgene3"){
    plot(density(vector), xlim=c(0, x), ylim=c(0, y), xlab="normalized tree distance", col=color)
  }
  else {
    lines(density(vector), col=color)
  }
}
 
pdf("bsppdist_mc_yule.pdf")
#### mc with yule model matrix
for (a in 1:300){
  #plotdist(bsgenes,a,5,"bsgenes","blue")
  plotdist(bsgene3,a, 5, "bsgene3","blue",3,15)
  #plotdist(bsgene2,a, 3, "bsgene2","red")
  plotdist(ppgene3,a, 5, "ppgene3","red",3,15)
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}
 
dev.off()

pdf("bsdist_mc_unif.pdf")
#### mc with yule model matrix
for (a in 1:300){
  #plotdist(bsgenes,a,6,"bsgenes","blue")
  #plotdist(bsgene2,a, 4, "ppgenes","red")
  plotdist(bsgene3,a, 6, "bsgene3","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("bsppdist_rc_yule.pdf")
#### mc with yule model matrix
for (a in 1:300){
 # plotdist(bsgenes,a,8,"bsgenes","blue")
  #plotdist(bsgene2,a, 6, "ppgenes","red")
  plotdist(bsgene3,a, 8, "bsgene3","blue",1.5,30)
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
  plotdist(ppgene3,a, 8, "ppgenes","red",1.5,30)
  
}

dev.off()

pdf("bspdist_rc_unif.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(bsgenes,a,9,"bsgenes","blue")
  plotdist(bsgene2,a, 7, "ppgenes","red")
  plotdist(bsgene3,a, 9, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("bsppdist_tt_yule.pdf")
#### mc with yule model matrix
for (a in 1:300){
  #plotdist(bsgenes,a,11,"bsgenes","blue")
  #plotdist(bsgene2,a, 9, "ppgenes","red")
  plotdist(bsgene3,a, 11, "bsgene3","blue",3,35)
  plotdist(ppgene3,a, 11, "ppgenes","red",3,35)
  
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("bsdist_tt_unif.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(bsgenes,a,12,"bsgenes","blue")
  plotdist(bsgene2,a, 10, "ppgenes","red")
  plotdist(bsgene3,a, 12, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

############for pp

pdf("ppdist_mc_yule.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(ppgenes,a,5,"bsgenes","blue")
  plotdist(ppgene3,a, 5, "ppgene3","black")
  plotdist(ppgene2,a, 3, "ppgene2","red")
  
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("ppdist_mc_unif.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(ppgenes,a,6,"bsgenes","blue")
  plotdist(ppgene2,a, 4, "ppgenes","red")
  plotdist(ppgene3,a, 6, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("ppdist_rc_yule.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(ppgenes,a,8,"bsgenes","blue")
  plotdist(ppgene2,a, 6, "ppgenes","red")
  plotdist(ppgene3,a, 8, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("ppdist_rc_unif.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(ppgenes,a,9,"bsgenes","blue")
  plotdist(ppgene2,a, 7, "ppgenes","red")
  plotdist(ppgene3,a, 9, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("ppdist_tt_yule.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(ppgenes,a,11,"bsgenes","blue")
  plotdist(ppgene2,a, 9, "ppgenes","red")
  plotdist(ppgene3,a, 11, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()

pdf("ppdist_tt_unif.pdf")
#### mc with yule model matrix
for (a in 1:300){
  plotdist(ppgenes,a,12,"bsgenes","blue")
  plotdist(ppgene2,a, 10, "ppgenes","red")
  plotdist(ppgene3,a, 12, "ppgenes","black")
  #plotdist(G3genes,a,5,"mlgenes","green")
  #plotdist(G5genes,a,5,"ppcongenes","black")
  #plotdist(G6genes,a,5,"G6genes","purple")
}

dev.off()
