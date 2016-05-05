library(ape)
library(diversitree)
entiretree<-read.nexus("ppmpchrpl_660.tre")
plot(ladderize(entiretree))
tips<- c("Coelegnathus_helena","Coelognathis_subradiatus", "Coelognathus_erythrurus", "Coelognathus_flavolineatus", "Coelognathus_radiatus", 
         "Coluber_constrictor", "Drymobius_margaritiferus", "Gonyosoma_boulengeri", "Gonyosoma_frenatum", "Gonyosoma_oxycephalum", "Gonyosoma_prasinum",
         "Gyalopion_canum", "Hapsidophrys_lineatus", "Hemorrhois_ravergieri", "Nerodia_sipedon", "Ptyas_mucosa", "Tantilla_coronata")
ingrouptree<-drop.tip(entiretree, tips)
plot(ladderize(ingrouptree))
write.tree(ladderize(ingrouptree), file="ingrouptree_660.txt")

 
tree <- read.tree("ingrouptree_660.txt")
state <- read.table("rat_geostates.txt", head=T, sep="\t")
head(state)
tipstate <- state[, 2]

names(tipstate) <- state[,1]
head(tipstate)
tree$tip.state <-tipstate

p<-starting.point.geosse(tree)
lik1<-make.geosse(tree, tree$tip.state, strict =F)
lik2<-constrain(lik1, sAB ~ 0)
ml1 <- find.mle(lik1, p)
p <- coef(ml1)
ml2 <- find.mle(lik2, p[argnames(lik2)])
round(rbind(full = coef(ml1), no.sAB = coef(ml2, TRUE)), 2)
anova_out<-anova(ml1, no.sAB = ml2)
anova_out
  p <- coef(ml1)
  prior <- make.prior.exponential(1/2)
  tmp <- mcmc(lik1, p, nsteps=100, prior=prior, w=1, print.every=0)
  w <- diff(sapply(tmp[2:8], quantile, c(0.025, 0.975)))
  mcmc1 <- mcmc(lik1, p, nsteps=20000, prior=prior, w=w)
  mcmc1diff <- with(mcmc1, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB, div.A=sA-xA, div.B=sB-xB))
  colMeans(mcmc1diff > 0)
  col1 <- c("red", "orange", "green","blue", "purple", "black", "gray")
  col2 <- col1[c(1,4,6)]
  mcmc1diff <- with(mcmc1, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB)) 
  par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
  profiles.plot(mcmc1[2:8], col.line=col1, xlab="", ylab="")
  legend("topright", argnames(lik1), col=col1, lty=1)
  profiles.plot(mcmc1diff, col.line=col2, xlab="", ylab="")
  legend("topright", colnames(mcmc1diff), col=col2, lty=1)
  title(xlab="rate", ylab="posterior probability density", outer=T, line=-1)



p <- coef(ml2)
prior <- make.prior.exponential(1/2)
tmp <- mcmc(lik2, p, nsteps=100, prior=prior, w=1, print.every=0)
w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))
mcmc2 <- mcmc(lik2, p, nsteps=20000, prior=prior, w=w)
mcmc2diff <- with(mcmc2, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB, div.A=sA-xA, div.B=sB-xB))
colMeans(mcmc2diff > 0)
col1 <- c("red", "orange", "blue", "purple", "black", "gray")
col2 <- col1[c(1,3,5)]
mcmc2diff <- with(mcmc2, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB)) 
par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
profiles.plot(mcmc2[2:7], col.line=col1, xlab="", ylab="")
legend("topright", argnames(lik2), col=col1, lty=1)
profiles.plot(mcmc2diff, col.line=col2, xlab="", ylab="")
legend("topright", colnames(mcmc2diff), col=col2, lty=1)
title(xlab="rate", ylab="posterior probability density", outer=T, line=-1)


par(mfrow=c(1,1))
posteriors<-cbind(mcmc1[, 1], mcmc1[, 9], mcmc2[, 8])

profiles.plot(posteriors[, 2:3], col.line=col1, xlab="posteriors", ylab="density")
legend("topleft", c("pp-full", "pp-no.sAB"), col=col1, lty=1)

par(mfrow=c(4,2))
for (i in 2:9){
  plot(mcmc1[, 1], mcmc1[,i])
}

par(mfrow=c(4,2))
for (j in 2:8){
  plot(mcmc2[, 1], mcmc2[,j])
}
save.image(file="ppmp660_geosse.RData")

