library(ape)
library(DDD)
library(laser)
library(TreePar)
library(TreeSim)
library(BAMMtools)
############RUN BAMM to test if their is rate shift on phylogeny
#tree<-read.tree("ratdppdiv2.tre")
tree<-read.tree("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal2.txt")
#tree<-read.tree("bsrfchrpltree.txt")

is.ultrametric(tree) # T

is.binary.tree(tree) # T
# Now to check min branch length:
min(tree$edge.length) # >0
#setBAMMpriors(tree)
#library(BAMMtools)
#tree1<-read.tree("ppmpchrpltree.txt")
#plot(ladderize(tree1))
tree<-ladderize(tree)
plot(tree)
edata <- getEventData(tree, eventdata = "/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/bamm_output/ppmp_cal2_event_data.txt", burnin=0.3, type="diversification")
#edata <- getEventData(tree, eventdata = "ppmp1_event_data3.txt", burnin=0.3)
summary(edata)
mcmcout <- read.csv("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/bamm_output/ppmp_cal2_mcmc_out.txt", header=T)
#mcmcout <- read.csv("ppmp1_mcmc_out3.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.3 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
shift_probs <- summary(edata)
plot.bammdata(edata, lwd=2, legend=T)
index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2)
addBAMMshifts(e2, cex=2)
#prior<-read.csv("ppmp2_prior_probs.txt")
prior<-read.csv("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/bamm_output/ppmp_cal2_prior_probs.txt")
priorshifts <- getBranchShiftPriors(tree, expectedNumberOfShifts = 0)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5,set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css, pal = "temperature")
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts =0)
par(mar= c(0,6,0,0))
q <- plot.bammdata(best, lwd = 2, tau = 0.001, breaksmethod = 'jenks', pal = "temperature", labels=T, cex = 0.8)
addBAMMshifts(best, cex=2.5)
addBAMMlegend(q, location = "topleft")
 
msc.set <- maximumShiftCredibility(edata, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)
 

 

#par(mfrow=c(1,1))
plotRateThroughTime(edata, ratetype="netdiv", ylim=c(0,0.8),intervalCol="gray80", avgCol="gray70", lwd = 8, opacity = 1, useMedian = F)
plotRateThroughTime(edata, node = 103, nodetype = "include",add = TRUE,ratetype="netdiv", ylim=c(0,0.8),intervalCol="orangered",lwd = 3, avgCol="orangered", opacity = 0.01, useMedian = F)
plotRateThroughTime(edata, node = 103, nodetype = "exclude",add = TRUE, ratetype="netdiv", ylim=c(0,0.8),intervalCol="steelblue",lwd = 3, avgCol="steelblue4", opacity = 0.01, useMedian = F)
legend("topright", legend = c("all ratsnakes", "shift-rate clade", "non shift-rate calde"), col = c("gray70", "orangered", "steelblue4"), lty = 1, lwd = 3)



par(mfrow=c(3,3))
plotRateThroughTime(edata, ratetype="speciation", ylim=c(0,0.8),intervalCol="gray50", avgCol="gray50")
plotRateThroughTime(edata, ratetype="extinction", ylim=c(0,0.8),intervalCol="gray50", avgCol="gray50")
plotRateThroughTime(edata, ratetype="netdiv", ylim=c(0,0.8),intervalCol="gray50",avgCol="gray50")
plotRateThroughTime(edata, ratetype="speciation",node = 103, nodetype = "include", ylim=c(0,0.8),intervalCol="orangered", avgCol="orangered")
plotRateThroughTime(edata, ratetype="extinction", node = 103, nodetype = "include", ylim=c(0,0.8),intervalCol="orangered", avgCol="orangered")
plotRateThroughTime(edata, ratetype="netdiv", node = 103, nodetype = "include", ylim=c(0,0.8),intervalCol="orangered",avgCol="orangered")
plotRateThroughTime(edata, ratetype="speciation",node = 103, nodetype = "exclude", ylim=c(0,0.8),intervalCol="steelblue4", avgCol="steelblue4")
plotRateThroughTime(edata, ratetype="extinction", node = 103, nodetype = "exclude", ylim=c(0,0.8),intervalCol="steelblue4", avgCol="steelblue4")
plotRateThroughTime(edata, ratetype="netdiv", node = 103, nodetype = "exclude", ylim=c(0,0.8),intervalCol="steelblue4",avgCol="steelblue4")


 



 
###########medusa
library(geiger)
res<-medusa(ladderize(tree))
print(names(res))
print(res$summary)
plot(res, cex = 0.5, label.offset = 1, edge.width = 2)


tiff("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/Figure6.tiff", width = 1000, height = 1600, pointsize = 40)
par(mar=c(2,0,0,0))  
plot(res, label.offset = 1, edge.width = 2)
dev.off()



#Appropriate  aicc-threshold for a tree of 70 tips is: 3.837418.

#Step 1: lnLik=-213.1524; aicc=428.334; model=yule
#Step 2: lnLik=-204.8555; aicc=415.8887; shift at node 103; model=yule; cut=stem; # shifts=1

#No significant increase in aicc score. Disregarding subsequent piecewise models.

#Model.ID Shift.Node Cut.At Model Ln.Lik.part         r epsilon     r.low    r.high
#1        1         71   node  yule    -110.216 0.0776653      NA 0.0534641 0.1083021
#2        2        103   stem  yule   -94.63941 0.2105960      NA 0.1498309 0.2859102
#
#



############################misfits models

out<-misfitssingle (tree, filename="ppmptreedivmodels.txt") 






############ Time-dependent


###############Density-dependent


#############Trait-dependent