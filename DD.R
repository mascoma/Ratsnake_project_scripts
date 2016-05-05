library(TreePar)
library(ape)
library(DDD)
tree<-read.tree("ppmpchrpltree.txt")
###DDD test
bd_0<-bd_ML(x,tdmodel=0) 
bd_1<-bd_ML(x,tdmodel=1)
bd_2<-bd_ML(x,tdmodel=2)
bd_3<-bd_ML(x,tdmodel=3)

dd_1<-dd_ML(x,ddmodel=1)
dd_2<-dd_ML(x,ddmodel=2)
dd_2_1<-dd_ML(x,ddmodel=2.1)
dd_2_2<-dd_ML(x,ddmodel=2.2)
dd_3<-dd_ML(x,ddmodel=3)
dd_4<-dd_ML(x,ddmodel=4)
dd_4.1<-dd_ML(x,ddmodel=4.1)

dd_SR_1<-dd_SR_ML(x,ddmodel=1)
dd_SR_2<-dd_SR_ML(x,ddmodel=2)
dd_SR_2_1<-dd_SR_ML(x,ddmodel=2.1)
dd_SR_2_2<-dd_SR_ML(x,ddmodel=2.2)
dd_SR_3<-dd_SR_ML(x,ddmodel=3)
dd_SR_4<-dd_SR_ML(x,ddmodel=4)
dd_SR_4.1<-dd_SR_ML(x,ddmodel=4.1)


x<-getx(tree)
bd.densdep.optim(x)

results: [1] 105
$par
[1] 0.4416676 0.1671704

$value
[1] 217.166

$counts
[1] 351

$convergence
[1] 0

$message
NULL

$hessian
NULL

$par
[1] 0.4416676 0.1671704

$value
[1] 217.166

$counts
[1] 351

$convergence
[1] 0

$message
NULL

$hessian
NULL

[[1]]
[[1]]$par
[1]   0.4416676   0.1671704 105.0000000

[[1]]$value
[1] 217.166

[[1]]$counts
[1] 351

[[1]]$convergence
[1] 0

[[1]]$message
NULL

[[1]]$hessian
NULL


[[2]]
[1] 0



#####treepar
timevec<-c(0,0.15,0.25)
lambdavec<-c(2.5,2,3)
muvec<-c(0.5,0.7,0.6)
x<-c(0.3,0.19,0.1)
x1<-c(x,max(x)*1.1)
x2<-c(x,max(x))
sampling<-0.4
grouptime<- rep(min(x)*0.95,length(x)+1)
group<- cbind(grouptime,grouptime*0+1)
group2 <- group
group2[1,2] <- 4
group2[2,2] <- 5
group2[3,2] <- 3
group3<-group
group3[2,2]<-10

### calculate likelihoods with root = 1

## Shifts in speciation / extinction rates (Stadler, PNAS 2011; Smrckova & Stadler, Manuscript 2014)
for (survival in c(0,1)) {
  print(LikShiftsPP(x,timevec,lambdavec,muvec,sampling,survival=survival))
  print(LikShifts(x,timevec,lambdavec,muvec,c(sampling,1,1),survival= survival))
  print(LikShifts(x,timevec,lambdavec,muvec,c(sampling,1,1),survival= survival,groups=group))
  print(LikShiftsSTT(par=c(lambdavec,muvec,timevec[-1]),x,x*0+1,sprob=c(0,0,0),
                     sampling=c(sampling,0,0),survival=survival,root=1))
  print(" ")
}

## Shifts in speciation / extinction rates with group sampling
for (survival in c(0,1)) {
  print(LikShifts(x,timevec,lambdavec,muvec,c(sampling,1,1),survival= survival,groups=group2))
  print(LikShifts(x,timevec,lambdavec,muvec,c(sampling,1,1),survival= survival,groups=group3))
  print(" ")
}

## Constant speciation and extinction rates 
# condition on age of tree x[1] and number of tips n
LikShiftsPP(x,timevec[1],lambdavec[1],muvec[1],sampling,n=1)
LikConstantn(lambdavec[1],muvec[1],sampling,x)
print(" ")
# condition on age of tree x[1]
for (survival in c(0,1)) {
  print(LikConstant(lambdavec[1],muvec[1],sampling,x,root=1,survival=survival))
  print(LikShiftsSTT(par=c(lambdavec[1],lambdavec[1],muvec[1],muvec[1],1),x,x*0+1,
                     sprob=c(0,0),sampling=c(sampling,1),survival=survival,root=1))
  print(LikShiftsPP(x,c(0),lambdavec[1],muvec[1],sampling,root=1,survival=survival))
  print(LikShifts(x,c(0),lambdavec[1],muvec[1],c(sampling),survival=survival))
  print(LikShifts(x,c(0),lambdavec[1],muvec[1],c(sampling),survival= survival,groups=group))
  print(LikShifts(x,c(0),c(lambdavec[1],lambdavec[1],lambdavec[1]),
                  c(muvec[1],muvec[1],muvec[1]),c(sampling,1,1),survival= survival))
  if (survival == 0 ) {
    print(LikDD(c(lambdavec[1],muvec[1], 200), 
                model=0 ,root=1, x=sort(x),sampling=sampling)[1])  }
  if (survival == 0 ) {
    print(LikDD(c(lambdavec[1],muvec[1], 300), 
                model=-1 ,root=1, x=sort(x),sampling=sampling)[1])  }
  print(" ")
}

## Diversity-dependent speciation rates
# condition on age of tree x[1], survival = 0
N<-10
pars <- matrix(c(N,lambdavec[1],muvec[1],0,sampling),nrow=1)
library(expoTree)
print(-runExpoTree(pars,sort(x2),rep(1,length(x2)),survival=0)+(length(x)-1)*log(2))
print(LikDD(c(lambdavec[1],muvec[1],N),x=sort(x),model=-1,root=1,sampling=sampling)[1])
print(" ")

### calculate likelihoods with root = 0

## Constant speciation and extinction rates 
# condition on age of tree x[1] and number of tips n
print(LikShiftsPP(x,timevec[1],lambdavec[1],muvec[1],sampling,root=0,n=1))
print(LikConstantn(lambdavec[1],muvec[1],sampling,x,root=0))
print(" ")
# condition on age of tree x[1]
for (survival in c(0,1)){
  print(LikShiftsPP(x,c(0),lambdavec[1],muvec[1],sampling,root=0,survival=survival))
  print(LikConstant(lambdavec[1],muvec[1],sampling,x,root=0,survival=survival))
  if (survival == 0 ) {print(LikDD(c(lambdavec[1],muvec[1], 200), 
                                   model=0 ,root=0, x=sort(x),sampling=sampling)[1])  }
  if (survival == 0 ) {print(LikDD(c(lambdavec[1],muvec[1], 300), 
                                   model=-1 ,root=0, x=sort(x),sampling=sampling)[1])  }
  print(" ")
}
## Diversity-dependent speciation rates
# condition on age of tree x[1], survival = 0
print(-runExpoTree(pars,sort(x),rep(1,length(x)),survival=0)+(length(x)-1)*log(2))
print(LikDD(c(lambdavec[1],muvec[1],N),x=sort(x),model=-1,root=0,sampling=sampling)[1])
