library(ape)
library(diversitree)
 

 
tree <- read.tree("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal2.txt")
state <- read.table("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/rat_geostates.txt", head=T, sep="\t")
head(state)
tipstate <- state[, 2]

names(tipstate) <- state[,1]
head(tipstate)
tree$tip.state <-tipstate

p<-starting.point.geosse(tree)
lik1<-make.geosse(tree, tree$tip.state, strict =F)
lik2<-constrain(lik1, sAB ~ 0)
lik3<-constrain(lik1, sA ~ sB, xA ~ xB, sAB~0)
ml1 <- find.mle(lik1, p)
p <- coef(ml1)
ml2 <- find.mle(lik2, p[argnames(lik2)])
ml3 <- find.mle(lik3, p[argnames(lik3)])
round(rbind(full = coef(ml1), no.sAB = coef(ml2, TRUE),eq.div = coef(ml3, TRUE)),3)
anova_out<-anova(ml1, no.sAB = ml2, eq.div = ml3)
anova_out2 <-anova(ml3,  diff.div = ml2)
anova_out
anova_out2
#Df   lnLik    AIC  ChiSq Pr(>|Chi|)   
#full    7 -212.57 439.15                     
#no.sAB  6 -214.31 440.63  3.480   0.062115 . 
#eq.div  4 -219.80 447.59 14.447   0.002355 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

  p <- coef(ml2)
  prior <- make.prior.exponential(1/2)
  set.seed(1)
  tmp <- mcmc(lik2, p, nsteps=100, prior=prior, w=1, print.every=0)
  w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))
  mcmc2 <- mcmc(lik2, p, nsteps=20000, prior=prior, w=w)
  mcmc2diff <- with(mcmc2, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB, div.A=sA-xA, div.B=sB-xB))
  colMeans(mcmc2diff > 0)
  col1 <- c("red", "orange", "green","blue", "purple", "black", "gray")
  col2 <- col1[c(1,4,6)]
  mcmc1diff <- with(mcmc1, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB)) 
  par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
  profiles.plot(mcmc1[2:8], col.line=col1, xlab="", ylab="")
  legend("topright", argnames(lik1), col=col1,pch = 15)
  profiles.plot(mcmc1diff, col.line=col2, xlab="", ylab="")
  legend("topright", colnames(mcmc1diff), col=col2, pch =15)
  title(xlab="rate", ylab="posterior probability density", outer=T, line=-1)



p <- coef(ml2)
prior <- make.prior.exponential(1/2)
tmp <- mcmc(lik2, p, nsteps=100, prior=prior, w=1, print.every=0)
w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))
mcmc2 <- mcmc(lik2, p, nsteps=20000, prior=prior, w=w)
burnstart <- floor(0.3 * nrow(mcmc2))
postburn <- mcmc2[burnstart:nrow(mcmc2), ]
library(coda)
effectiveSize(postburn$p)
 
effectiveSize(postburn$dA)

mcmc2diff <- with(postburn, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB, div.A=sA-xA, div.B=sB-xB))
#colMeans(mcmc2diff > 0)
col1 <- c("red", "orange", "blue", "purple", "gray20", "gray70")
col2 <- c("deeppink4","forestgreen","gray10", "darkorange4")
mcmc2diff <- with(postburn, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB, div.diff=((sA-xA)-(sB-xB)))) 
#par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
profiles.plot(postburn[2:7], col.line=col1, xlab="rate", ylab="posterior probability density")
legend("topright", argnames(lik2), col=col1, lty = 1, lwd = 2)
profiles.plot(mcmc2diff, col.line=col2, xlab="", ylab="")
legend("topright", c("sA-sB", "xA-xB", "dA-dB", "divA-divB"), col=col2, lty=1, lwd=2)
title(xlab="rate", ylab="posterior probability density", outer=T, line=-1)
s.diff=postburn$sA-postburn$sB
x.diff=postburn$xA-postburn$xB
d.diff=postburn$dA-postburn$dB
div.A=postburn$sA-postburn$xA
div.B=postburn$sB-postburn$xB
div.diff = div.A - div.B
t.test(s.diff)
t.test(x.diff) 
t.test(d.diff)
t.test(div.diff)


par(mfrow=c(4,2))
for (i in 2:9){
  plot(mcmc1[, 1], mcmc1[,i])
}

par(mfrow=c(4,2))
for (j in 2:8){
  plot(mcmc2[, 1], mcmc2[,j])
}
save.image(file="/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/geosse.RData")

