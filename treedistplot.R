##to compare gene tree distribution of different groups 
### using 300 replicates 
### graph 1 will be mc matrix, yule average, 6 groups
### graph 2 will be mc matrix, unifor average, 6 groups
### graph 3 will be RF matrix, yule average, 6 groups
### graph 4 will be RF matrix, unifor average, 6 groups
### graph 5 will be tt matrix, yule average, 6 groups
### group 6 will be tt matrix, unifor average, 6 groups
 G1file1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G1bs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
 G2file1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G2bs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
 G3file1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G3bs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
 G5file1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G5bs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)
 G6file1 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G6bs1000outforgraph", pattern="*.txt", full.names=T, recursive=F)

readtable2<-function (file) {read.table(file,header=T)}

 G1genes<-lapply(G1file1, function(x) readtable2(x))
 G2genes<-lapply(G2file1, function(x) readtable2(x))
 G3genes<-lapply(G3file1, function(x) readtable2(x))
 G5genes<-lapply(G5file1, function(x) readtable2(x))
 G6genes<-lapply(G6file1, function(x) readtable2(x))
 
 
 
 ###################################
 g1<-G1genes[[10]]
 g2<-G2genes[[10]]
 g3<-G3genes[[10]]
 g5<-G5genes[[10]]
 g6<-G6genes[[10]]
 bs<-bsgenes[[10]]
 
 bsmc<-bs[[5]]
 plot(density(log(bsmc)))
 sbsmc<-sqrt(bsmc)
 
 g2mc<-g2[[5]]
 plot(density(sqrt(g2mc)))
 qqnorm(sqrt(g2mc))
 sg2mc<-sqrt(g2mc)
 
 
 g3mc<-g3[[5]]
 plot(density(sqrt(g3mc)))
 qqnorm(g3mc)
 sg3mc<-sqrt(g3mc)
 
 
 g5mc<-g5[[5]]
 plot(density(sqrt(g5mc)))
 qqnorm(sqrt(g5mc))
 sg5mc<-sqrt(g5mc)
 
 
 g6mc<-g6[[5]]
 plot(density(sqrt(g6mc)))
 qqnorm(sqrt(g6mc)) 
 sg6mc<-sqrt(g6mc)
 
 
 
 plot(density(sg2mc), ylim=c(0,6))
 lines(density(sg3mc))
 lines(density(sg5mc))
 lines(density(sg6mc))
 line(density(bsmc))
 t.test(lbsmc,lppmc)
 
 fac<-rep(c("g2","g3","g5","g6","bs"),each=46056)
 v1<-c(sg2mc,sg3mc,sg5mc,sg6mc,bsmc)
 mcyule<-data.frame(v1,fac)
 head(mcyule)
 
 result<-aov(v1~fac,data=mcyule)
 summary(result)
 TukeyHSD(result,conf.level=0.95)
 pairwise.t.test(v1, fac, p.adjust="bonferroni")
 plot(v1~fac, data=mcyule)
 
 bsrc<-bs[[8]]
 plot(density(sqrt(bsrc)))
 qqnorm(sqrt(bsrc))
 sbsrc<-sqrt(bsrc)
 
 g2rc<-g2[[8]]
 plot(density(sqrt(g2rc)))
 qqnorm(sqrt(g2rc))
 sg2rc<-sqrt(g2rc)
 
 
 g3rc<-g3[[8]]
 plot(density(sqrt(g3rc)))
 qqnorm(sqrt(g3rc))
 sg3rc<-sqrt(g3rc)
 
 
 g5rc<-g5[[8]]
 plot(density(sqrt(g5rc)))
 qqnorm(sqrt(g5rc))
 sg5rc<-sqrt(g5rc)
 
 
 g6rc<-g6[[8]]
 plot(density(sqrt(g6rc)))
 qqnorm(sqrt(g6rc)) 
 sg6rc<-sqrt(g6rc)
 
 
 
 plot(density(sg2rc), ylim=c(0,35))
 lines(density(sg3rc))
 lines(density(sg5rc))
 lines(density(sg6rc))
 lines(density(sbsrc))
# t.test(lbsmc,lppmc)
 
 fac<-rep(c("g2","g3","g5","g6","bs"),each=46056)
 v1<-c(sg2rc,sg3rc,sg5rc,sg6rc,sbsrc)
 rcyule<-data.frame(v1,fac)
 head(rcyule)
 
 result<-aov(v1~fac,data=rcyule)
 summary(result)
 TukeyHSD(result,conf.level=0.95)
 pairwise.t.test(v1, fac, p.adjust="bonferroni")
 
 plot(v1~fac, data=rcyule)
  
 
 bstt<-bs[[11]]
 plot(density(log(1/bstt)))
 lbstt<-log(1/bstt)
 
 g2tt<-g2[[11]]
 plot(density(log(1/(g2tt))))
 qqnorm(log(1/(g2tt)))
 lg2tt<-log(1/g2tt)
 
 
 g3tt<-g3[[11]]
 plot(density(log(1/(g3tt))))
 qqnorm(log(1/g3tt))
 lg3tt<-log(1/g3tt)
 
 
 g5tt<-g5[[11]]
 plot(density(log(1/(g5tt))))
 qqnorm(log(1/(g5tt)))
 lg5tt<-log(1/(g5tt))
 
 
 g6tt<-g6[[11]]
 plot(density(log(1/(g6tt))))
 qqnorm(log(1/(g6tt))) 
 lg6tt<-log(1/(g6tt))
 
 
 
 plot(density(lg2tt), ylim=c(0,3))
 lines(density(lg3tt))
 lines(density(lg5tt))
 lines(density(lg6tt))
 lines(density(lbstt))
 # t.test(lbsmc,lppmc)
 
 fac<-rep(c("g2","g3","g5","g6","bs"),each=46056)
 v1<-c(lg2tt,lg3tt,lg5tt,lg6tt,lbstt)
 ttyule<-data.frame(v1,fac)
 head(ttyule)
 
 result<-aov(v1~fac,data=ttyule)
 summary(result)
 TukeyHSD(result,conf.level=0.95)
 pairwise.t.test(v1, fac, p.adjust="bonferroni")
 
 plot(v1~fac, data=ttyule)
 
 
 
 
 
 
 #wmatrix<- G1genes[[1]]
 #vector<-matrix[,5]
# plot(density(vector), xlim=c(0, 1.5), ylim=c(0, 3), xlab="normalized tree distance", col="blue")
 #file<-G1genes
 #names(file)
 plotdist<-function(files, a, b, name, color, x,y){
  matrix<-files[[a]]
  vector<-matrix[,b] 
  if (a==1 && name == "G2genes"){
   plot(density(vector), xlim=c(0, x), ylim=c(0,y), xlab="normalized tree distance", col=color)
  }
  else {
    lines(density(vector), col=color)
  }
}

 
 
 pdf("bssubcladesgenedist_mc_yule.pdf")
 #### gene distribution mc with unifor model matrix 
 for (a in 1:300){
   #plotdist(G1genes,a,5,"G1genes","blue")
   plotdist(G2genes,a, 5, "G2genes","red",2,4)
   plotdist(G3genes,a,5,"G3genes","blue",2,4)
   plotdist(G5genes,a,5,"G5genes","black",2,4)
   plotdist(G6genes,a,5,"G6genes","purple",2,4)
   plotdist(bsgenes,a,5,"bsgenes","green",2,4)
 }
 dev.off()
 
 pdf("genedist_mc_unif.pdf")
#### gene distribution mc with unifor model matrix 
 for (a in 1:300){
   #plotdist(G1genes,a,6,"G1genes","blue")
   plotdist(G2genes,a, 6, "G2genes","red")
   plotdist(G3genes,a,6,"G3genes","blue")
   plotdist(G5genes,a,6,"G5genes","black")
   plotdist(G6genes,a,6,"G6genes","purple")
   
 }
 dev.off()
 
 pdf("bssubcladesgenedist_rc_yule.pdf")
### rc with yule model 
 for (a in 1:300){
  #plotdist(G1genes,a,8,"G1genes","blue", 20)
   plotdist(G2genes,a, 8, "G2genes","red",1.5,25)
   plotdist(G3genes,a,8,"G3genes","blue",1.5,25)
   plotdist(G5genes,a,8,"G5genes","black",1.5,25)
   plotdist(G6genes,a,8,"G6genes","purple",1.5,25)
   plotdist(bsgenes,a,8,"bsgenes","green",1.5,25)
 }
 dev.off()
 
 pdf("bssubcladesgenedist_rc_unif.pdf")
 
### rc with unif model 
 for (a in 1:300){
   plotdist(G1genes,a,9,"G1genes","blue")
   plotdist(G2genes,a, 9, "G2genes","red")
   plotdist(G3genes,a,9,"G3genes","green")
   plotdist(G5genes,a,9,"G5genes","black")
   plotdist(G6genes,a,9,"G6genes","purple")
 }
 dev.off()
 
 pdf("bssubcladesgenedist_tt_yule.pdf")
### ttiwith yule model
 
 for (a in 1:300){
   #plotdist(G1genes,a,11,"G1genes","blue")
   plotdist(G2genes,a, 11, "G2genes","red",1.5,4)
   plotdist(G3genes,a,11,"G3genes","blue",1.5,4)
   plotdist(G5genes,a,11,"G5genes","black",1.5,4)
   plotdist(G6genes,a,11,"G6genes","purple",1.5,4)
   plotdist(bsgenes,a,11,"bsgenes","green",1.5,4)
 }
 dev.off()
 
 pdf("genedist_tt_unif.pdf")
 #####tt with unifr model
 for (a in 1:300){
   plotdist(G1genes,a,12,"G1genes","blue")
   plotdist(G2genes,a, 12, "G2genes","red")
   plotdist(G3genes,a,12,"G3genes","green")
   plotdist(G5genes,a,12,"G5genes","black")
   plotdist(G6genes,a,12,"G6genes","purple")
 }
 dev.off()
 
 
###########################################################3
 ####################################################
 #################################################
 ###############################################
 ################ compare gene trees vs specids trees
 
 G1file2 <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G1gtvssptoutput", pattern="*.txt", full.names=T, recursive=F)
 G2mp <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G2gtvssptoutput/mptreedist", pattern="*.txt", full.names=T, recursive=F)
 G2RF <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G2gtvssptoutput/RFtreedist", pattern="*.txt", full.names=T, recursive=F)
 
 G3mp <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G3gtvssptoutput/mptreedist", pattern="*.txt", full.names=T, recursive=F)
 G3RF <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G3gtvssptoutput/RFtreedist", pattern="*.txt", full.names=T, recursive=F)
 
 G5mp <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G5gtvssptoutput/mptreedist", pattern="*.txt", full.names=T, recursive=F)
 G5RF <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G5gtvssptoutput/RFtreedist", pattern="*.txt", full.names=T, recursive=F)
 
 
 G6mp <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G6gtvssptoutput/mptreedist", pattern="*.txt", full.names=T, recursive=F)
 G6RF <- list.files(path="G:/work/treedistribution/treecmp/sumresults/G6gtvssptoutput/RFtreedist", pattern="*.txt", full.names=T, recursive=F)
 
 readtable2<-function (file) {read.table(file,header=T)}
 
 G1mp<-lapply(G1file2, function(x) readtable2(x))
 G2mpt<-lapply(G2mp, function(x) readtable2(x))
 G2RFt<-lapply(G2RF, function(x) readtable2(x))
 G3mpt<-lapply(G3mp, function(x) readtable2(x))
 G3RFt<-lapply(G3RF, function(x) readtable2(x))
 
 G5mpt<-lapply(G5mp, function(x) readtable2(x))
 G5RFt<-lapply(G5RF, function(x) readtable2(x))
 
 G6mpt<-lapply(G6mp, function(x) readtable2(x))
 G6RFt<-lapply(G6RF, function(x) readtable2(x))
 
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
 pdf("g2gtvsspt_mc_yule.pdf")
 for (a in 1:300){
   plotdist(G2mpt, a, 3, "G2mpt","G2mpt","blue",2,4)
   plotdist(G2RFt, a, 3, "G2mpt", "G2RFt", "red",2,4)
   plotdist(G2genes, a, 5, "G2mpt","G2genes","black",2,4)
 }
 
 dev.off()
 
 
 pdf("g3gtvsspt_mc_yule.pdf")
 for (a in 1:300){
   plotdist(G3mpt, a, 3, "G3mpt","G3mpt","blue",2,4)
  # plotdist(G3RFt, a, 3, "G3mpt", "G3RFt", "red",2,4)
   plotdist(G3genes, a, 5, "G3mpt","G3genes","black",2,4)
 }
 
 dev.off()
 
 pdf("g5gtvsspt_mc_yule.pdf")
 for (a in 1:300){
   plotdist(G5mpt, a, 3, "G5mpt","G5mpt","blue",2,4)
  # plotdist(G5RFt, a, 3, "G5mpt", "G5RFt", "red",2,4)
   plotdist(G5genes, a, 5, "G5mpt","G5genes","black",2,4)
 }
 
 dev.off()
 
 
 pdf("g6gtvsspt_mc_yule.pdf")
 for (a in 1:300){
   plotdist(G6mpt, a, 3, "G6mpt","G6mpt","blue",2,4)
   #plotdist(G6RFt, a, 3, "G6mpt", "G6RFt", "red",2,4)
   plotdist(G6genes, a, 5, "G6mpt","G6genes","black",2,4)
 }
 
 dev.off()
 
 ##################### RC
 
 pdf("g2gtvsspt_rc_yule.pdf")
 for (a in 1:300){
   plotdist(G2mpt, a, 6, "G2mpt","G2mpt","blue",1.5,20)
   plotdist(G2RFt, a, 6, "G2mpt", "G2RFt", "red",1.5,20)
   plotdist(G2genes, a, 8, "G2mpt","G2genes","black",1.5,20)
 }
 
 dev.off()
 
 
 pdf("g3gtvsspt_rc_yule.pdf")
 for (a in 1:300){
   plotdist(G3mpt, a, 6, "G3mpt","G3mpt","blue",1.5,25)
   # plotdist(G3RFt, a, 3, "G3mpt", "G3RFt", "red",2,4)
   plotdist(G3genes, a, 8, "G3mpt","G3genes","black",1.5,25)
 }
 
 dev.off()
 
 pdf("g5gtvsspt_rc_yule.pdf")
 for (a in 1:300){
   plotdist(G5mpt, a, 6, "G5mpt","G5mpt","blue",1.5,15)
   # plotdist(G5RFt, a, 3, "G5mpt", "G5RFt", "red",2,4)
   plotdist(G5genes, a, 8, "G5mpt","G5genes","black",1.5,15)
 }
 
 dev.off()
 
 
 pdf("g6gtvsspt_rc_yule.pdf")
 for (a in 1:300){
   plotdist(G6mpt, a, 6, "G6mpt","G6mpt","blue",1.5,25)
   #plotdist(G6RFt, a, 3, "G6mpt", "G6RFt", "red",2,4)
   plotdist(G6genes, a, 8, "G6mpt","G6genes","black",1.5,25)
 }
 
 dev.off()
 
 
##########################tt
 
 
 
 pdf("g2gtvsspt_tt_yule.pdf")
 for (a in 1:300){
   plotdist(G2mpt, a, 9, "G2mpt","G2mpt","blue",1.5,4)
   plotdist(G2RFt, a, 9, "G2mpt", "G2RFt", "red",1.5,4)
   plotdist(G2genes, a, 11, "G2mpt","G2genes","black",1.5,4)
 }
 
 dev.off()
 
 
 pdf("g3gtvsspt_tt_yule.pdf")
 for (a in 1:300){
   plotdist(G3mpt, a, 9, "G3mpt","G3mpt","blue",1.5,4)
   # plotdist(G3RFt, a, 3, "G3mpt", "G3RFt", "red",2,4)
   plotdist(G3genes, a, 11, "G3mpt","G3genes","black",1.5,4)
 }
 
 dev.off()
 
 pdf("g5gtvsspt_tt_yule.pdf")
 for (a in 1:300){
   plotdist(G5mpt, a, 9, "G5mpt","G5mpt","blue",1.5,4)
   # plotdist(G5RFt, a, 3, "G5mpt", "G5RFt", "red",2,4)
   plotdist(G5genes, a, 11, "G5mpt","G5genes","black",1.5,4)
 }
 
 dev.off()
 
 
 pdf("g6gtvsspt_tt_yule.pdf")
 for (a in 1:300){
   plotdist(G6mpt, a, 9, "G6mpt","G6mpt","blue",1.5,4)
   #plotdist(G6RFt, a, 3, "G6mpt", "G6RFt", "red",2,4)
   plotdist(G6genes, a,11, "G6mpt","G6genes","black",1.5,4)
 }
 
 dev.off()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 pdf("gtspt_mc_yule.pdf")
 #### mc with yule model matrix
 
 head(G1gene2[[1]])
 for (a in 1:600){
   plotdist(G1gene2,a,3,"G1gene2","blue")
   plotdist(G2gene2,a, 3, "G2gene2","red")
   plotdist(G3gene2,a,3,"G3gene2","green")
   plotdist(G5gene2,a,3,"G5gene2","black")
   plotdist(G6gene2,a,3,"G6gene2","purple")
 }
dev.off()
 
 pdf("gtspt_mc_unif.pdf")
 #### mc with unifor model matrix 
 for (a in 1:600){
   plotdist(G1gene2,a,4,"G1gene2","blue")
   plotdist(G2gene2,a, 4, "G2gene2","red")
   plotdist(G3gene2,a,4,"G3gene2","green")
   plotdist(G5gene2,a,4,"G5gene2","black")
   plotdist(G6gene2,a,4,"G6gene2","purple")
 }
 dev.off()
 
 pdf("gtspt_rc_yule.pdf")
 ### rc with yule model 
 for (a in 1:600){
   plotdist(G1gene2,a,6,"G1gene2","blue")
   plotdist(G2gene2,a, 6, "G2gene2","red")
   plotdist(G3gene2,a,6,"G3gene2","green")
   plotdist(G5gene2,a,6,"G5gene2","black")
   plotdist(G6gene2,a,6,"G6gene2","purple")
 }
 dev.off()
 
 pdf("gtsptrc_unif.pdf")
 
 ### rc with unif model 
 for (a in 1:600){
   plotdist(G1gene2,a,7,"G1gene2","blue")
   plotdist(G2gene2,a, 7, "G2gene2","red")
   plotdist(G3gene2,a,7,"G3gene2","green")
   plotdist(G5gene2,a,7,"G5gene2","black")
   plotdist(G6gene2,a,7,"G6gene2","purple")
 }
 dev.off()
 
 pdf("gtspt_tt_yule.pdf")
 ### ttiwith yule model
 
 for (a in 1:600){
   plotdist(G1gene2,a,9,"G1gene2","blue")
   plotdist(G2gene2,a, 9, "G2gene2","red")
   plotdist(G3gene2,a,9,"G3gene2","green")
   plotdist(G5gene2,a,9,"G5gene2","black")
   plotdist(G6gene2,a,9,"G6gene2","purple")
 }
 dev.off()
 
 pdf("gtspt_tt_unif.pdf")
 #####tt with unifr model
 for (a in 1:600){
   plotdist(G1gene2,a,10,"G1gene2","blue")
   plotdist(G2gene2,a, 10, "G2gene2","red")
   plotdist(G3gene2,a,10,"G3gene2","green")
   plotdist(G5gene2,a,10,"G5gene2","black")
   plotdist(G6gene2,a,10,"G6gene2","purple")
 }
 dev.off()