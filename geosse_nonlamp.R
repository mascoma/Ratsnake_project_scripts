library(ape)
library(diversitree)

treedir <- "/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal2.txt"
statedir <- "/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/rat_geostates.txt"
outputdir <- "/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/"

tree.all <- read.tree(treedir)
state <- read.table(statedir, head=T, sep="\t")
head(state)
tipstate <- state[, 2]


names(tipstate) <- state[, 1]
head(tipstate)
tree.all$tip.state <- tipstate

tip.lamp <- tree.all$tip.label[49:70]
tip.nonlamp <- tree.all$tip.label[1:48]
tree.lamp <- drop.tip(tree.all, tip.nonlamp)
tree.nonlamp <- drop.tip(tree.all, tip.lamp)
tree.lamp$tip.state
tree.nonlamp$tip.state <- tree.all$tip.state[1:48]

tree <- tree.nonlamp
tree
tree$tip.state

p <- starting.point.geosse(tree)
lik1 <- make.geosse(tree, tree$tip.state, strict =F)
lik2 <- constrain(lik1, sAB ~ 0)
lik3 <- constrain(lik1, sA ~ sB, xA ~ xB, sAB~0)
ml1 <- find.mle(lik1, p)
p <- coef(ml1)
ml2 <- find.mle(lik2, p[argnames(lik2)])
ml3 <- find.mle(lik3, p[argnames(lik3)])
round(rbind(full = coef(ml1), no.sAB = coef(ml2, TRUE),eq.div = coef(ml3, TRUE)), 3)
anova_out <-anova(ml1, no.sAB = ml2, eq.div = ml3)
anova_out2 <- anova(ml3,  diff.div = ml2)
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
tmp <- mcmc(lik = lik2, x.init = p, nsteps = 100, prior = prior, w = 1, print.every = 0)
w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))
mcmc2 <- mcmc(lik2, p, nsteps = 20000, prior = prior, w = w)
burnstart <- floor(0.9 * nrow(mcmc2))
postburn <- mcmc2[burnstart:nrow(mcmc2), ]
library(coda)
ees.p <- effectiveSize(postburn$p)
ees.sA <- effectiveSize(postburn$sA)

mcmc2diff <- with(postburn, data.frame(s.diff = sA - sB, x.diff = xA - xB,
                  d.diff = dA - dB, div.diff= ((sA-xA) - (sB-xB))))
#colMeans(mcmc2diff > 0)
col1 <- c("red", "orange", "blue", "purple", "gray20", "gray70")
col2 <- c("deeppink4","forestgreen","gray10", "darkorange4")
#par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
image1 <- paste(outputdir, "nonlamp_allpar.png", sep = "")
png(image1)
profiles.plot(postburn[2:7], col.line=col1, xlab="rate", 
              ylab="posterior probability density")
legend("topright", argnames(lik2), col=col1, lty = 1, lwd = 2)
dev.off()
image2 <- paste(outputdir,"nonlamp_diff.png", sep = "")
png(image2)
profiles.plot(mcmc2diff, col.line=col2, xlab="", 
              ylab="")
legend("topright", c("sA-sB", "xA-xB", "dA-dB", "divA-divB"), col=col2, lty=1, 
       lwd=2)
title(xlab="rate", ylab="posterior probability density", outer=T, line=-1)
dev.off()

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
outputdata <- paste(outputdir,"geosse_nonlamp.RData", sep = "")
save.image(file = outputdata)


