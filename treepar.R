####treepar test

library(TreePar)
library(ape)

tree<-read.tree("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal2.txt")
x<-getx(tree)
bddens<-bd.densdep.optim(x)   ## density dependent model
bdyule<-bd.densdep.optim(x,Yule=T) ## density dependent with no extinction
bdshifts_1<-bd.shifts.optim(x,grid=1, sampling=c(1,1), start=0.5, end=35) ## birth-death model with one shift
bdshiftsyule_1<-bd.shifts.optim(x,grid=1, sampling=c(1,1), start=0.5, end=35, yule=T) ## birthdeath model with one shift, no extinction
bdshifts_2<-bd.shifts.optim(x,grid=1, sampling=c(1,1,1), start=0.5, end=35) ## birthdeath model with two shifts
bdshiftsyule_2<-bd.shifts.optim(x,grid=1, sampling=c(1,1,1), start=0.5, end=35, yule=T) #### birthdeath model with two shifts, no extinction

save.image(file= "ppmptreepar.RData")
quit()
load("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/treepar/ppmp_cal2_treepar.RData")
###########results   
## to determine best model using Akaike weights
## calculate AIC first
AIC_bddens<-2*length(bddens[[1]]$par)+2*bddens[[1]]$value
AIC_bdyule<-2*length(bdyule[[1]]$par)+2*bdyule[[1]]$value
AIC_bdshifts_0 <-2*2+2*2.131524e+02
AIC_bdshifts_1<- 2*5+2*2.102069e+02
AIC_bdshifts_2<- 2*8+2*208.99785224 
AIC_bdshiftsyule_0 <- 2*1+2*213.1523890
AIC_bdshiftsyule_1 <- 2*3+2*210.21026472
AIC_bdshiftsyule_2 <- 2*5+2*209.09996239

AICs<-c(AIC_bddens,AIC_bdyule,AIC_bdshifts_0,AIC_bdshifts_1,AIC_bdshifts_2,AIC_bdshiftsyule_0,AIC_bdshiftsyule_1,AIC_bdshiftsyule_2)
library(qpcR) 
akaike.weights(AICs)
#summary of the results Feb 18, 2016 
#$deltaAIC
#[1] 2.000000 0.000000 5.541452 5.650452 9.232357 3.541430 1.657182 3.436577

#$rel.LL
#[1] 0.367879441 1.000000000 0.062616517 0.059295243 0.009890521 0.170211211 0.436664150 0.179372867

#$weights
#[1] 0.160932071 0.437458724 0.027392141 0.025939221 0.004326695 0.074460379 0.191022542 0.078468226

 
 