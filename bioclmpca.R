
bioclm<-read.table("/Users/Xin/Documents/Ratsnake_project/input/bioclm_all.txt", header=T)
head(bioclm)
bc<-cbind(bioclm[,2:3],bioclm[,5:7],bioclm[,9:20])
head(bc)
fit<- princomp(bc,cor=T)
loadings(fit)

 summary(fit)
 head(fit$score)
bioclmpc_all <- cbind(as.character(bioclm[,1]), fit$score)
write.table(bioclmpc_all,file="/Users/Xin/Documents/Ratsnake_project/output/bioclmpca1.txt",sep="\t")

library(plotrix)

biopcs<-read.table("/Users/Xin/Documents/Ratsnake_project/output/bioclmpca1.txt",sep="\t",header=T)
head(biopcs)
mean1<-function(x){mean(x,na.rm=T)}
colnames <- vector("character", length = 35)
meanse<-matrix(,79,35)
meanse[,1]<-levels(biopcs[,1]) 
colnames[1] <- "species"

for(i in seq(2,35,by=2)){
  temp<-tapply(biopcs[,((i/2)+1)], biopcs[,1], FUN= mean1 )
  meanse[,i]<-temp
  colname = paste("mean", i/2, sep = "_")
  colnames[i]<- colname
}

for(j in seq(3, 35, by=2)){
  temp<-tapply(biopcs[, (((j-1)/2)+1)], biopcs[,1], FUN= std.error )
  meanse[,j]<-temp
  colname = paste("se", (j-1)/2, sep = "_")
  colnames[j]<- colname
}



meanse<-as.data.frame(meanse)
names(meanse)<-colnames

head(meanse)
write.table(meanse,file="/Users/Xin/Documents/Ratsnake_project/output/rat79biopcmeanse.txt",sep="\t")

