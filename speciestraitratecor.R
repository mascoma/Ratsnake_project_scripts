library(ape)
library(geiger)
tree<-read.tree("ppmptree66.tre")
trait<-read.table("tree66seprate.txt",header=T)
head(trait)
spe<-trait[,2]

trait1<-trait[,14]
#ed<-trait[,4]
#gape<-trait[,5]
#biopc4<-trait[,6]
names(trait1)<-names(spe)<-trait[,1]
treedata<-treedata(tree,spe)
phy<-treedata$phy
spe<-treedata$data[phy$tip.label,] 
Contrastspe<- pic(log(spe), phy)

treedata<-treedata(tree,trait1)  
phy<-treedata$phy
trait1<-treedata$data[phy$tip.label,] 
Contrasttrait1<- pic(log(trait1), phy)
spetrait1 <- lm(Contrastspe~Contrasttrait1)
summary.lm(spetrait1)
plot(Contrastspe, Contrasttrait1)
abline(spetrait1)
#bvf 0.43979      0.3689 
#ed 1.11e-10      0.359
#gape 5.57e-06 ***  0.2084 
#bp 6.61e-09 *** 0.3897
#srw1 0.0405 * 0.8965 
#srw2 5.51e-08 ***  0.6188 
#srw3 0.57876   0.06397
#svl 1.35e-05 *** 0.05471
#tl 3.48e-08 *** 0.07403
#drw1 0.00132 **  0.3092 
#drw2 0.40446        0.01797 *
#drw3 4.23e-06 *** 0.1404



tree<-read.tree("ppmptree69.tre")
trait<-read.table("tree69seprate.txt",header=T)

tree<-read.tree("ppmptree66.tre")
trait<-read.table("tree66seprate.txt",header=T)
head(trait)
species<-trait[,1]
data<-trait[,2:6]
rownames(data)<-trait[,1]
treedata<-treedata(tree,data)
phy<-treedata$phy
dat<-treedata$data[phy$tip.label,] 
 
species<-as.character(species)
dat<-cbind(species,dat)
dat<-as.data.frame(dat)
ratsnake<-comparative.data(phy,dat,species)
mod<-pgls(log(as.numeric(spe))~log(as.numeric(biopc4)),data=ratsnake,lambda='ML',k="ML",delta="ML") # lambda=1
summary(mod)
plot(mod)

spe<-trait[,2]
biopc1<-trait[,3]
biopc2<-trait[,4]
biopc3<-trait[,5]
biopc4<-trait[,6]
names(biopc4)<-names(biopc3)<-names(biopc2)<-names(biopc3)<-names(spe)<-trait[,1]
treedata<-treedata(tree,spe)
phy<-treedata$phy
spe<-treedata$data[phy$tip.label,] 
Contrastspe<- pic(log(spe), phy)



treedata<-treedata(tree,biopc1) #1.51e-05 ***  0.0048010 ** 
phy<-treedata$phy
biopc1<-treedata$data[phy$tip.label,] 
Contrastbiopc1<- pic(log(biopc1), phy)
spebio1 <- lm(Contrastspe~Contrastbiopc1)
summary.lm(spebio1)
plot(Contrastspe, Contrastbiopc1)
abline(spebio1)

treedata<-treedata(tree,biopc2) # 2.4e-11 *** 0.9174
phy<-treedata$phy
biopc2<-treedata$data[phy$tip.label,] 
Contrastbiopc2<- pic(log(biopc2), phy)
spebio2 <- lm(Contrastspe~Contrastbiopc2)
summary.lm(spebio2)
plot(Contrastspe, Contrastbiopc2)
abline(spebio2)

treedata<-treedata(tree,biopc3) #3.42e-10 ***  0.2131
phy<-treedata$phy
biopc3<-treedata$data[phy$tip.label,] 
Contrastbiopc3<- pic(log(biopc3), phy)
spebio3 <- lm(Contrastspe~Contrastbiopc3)
summary.lm(spebio3)
plot(Contrastspe, Contrastbiopc3)
abline(spebio3)

treedata<-treedata(tree,biopc4) #8.22e-08 ***   0.114  
phy<-treedata$phy
biopc4<-treedata$data[phy$tip.label,] 
Contrastbiopc4<- pic(log(biopc4), phy)
spebio4 <- lm(Contrastspe~Contrastbiopc4)
summary.lm(spebio4)
plot(Contrastspe, Contrastbiopc4)
abline(spebio4)