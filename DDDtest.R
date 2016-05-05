library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ingrouptree.txt")
x<-getx(tree)
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
dd_4.2<-dd_ML(x,ddmodel=4.2)
dd_5<-dd_ML(x,ddmodel=5)

dd_SR_1<-dd_SR_ML(x,ddmodel=1)
dd_SR_2<-dd_SR_ML(x,ddmodel=2)
dd_SR_2_1<-dd_SR_ML(x,ddmodel=2.1)
dd_SR_2_2<-dd_SR_ML(x,ddmodel=2.2)
dd_SR_3<-dd_SR_ML(x,ddmodel=3)
dd_SR_4<-dd_SR_ML(x,ddmodel=4)
dd_SR_4.1<-dd_SR_ML(x,ddmodel=4.1)
dd_SR_4.2<-dd_SR_ML(x,ddmodel=4.2)
save.image(file= "ppmpDDD.RData")
quit()
AIC_bd_0<-2*bd_0$df-2*bd_0$loglik
AIC_bd_1<-2*bd_1$df-2*bd_1$loglik
AIC_bd_2<-2*bd_2$df-2*bd_2$loglik
AIC_bd_3<-2*bd_3$df-2*bd_3$loglik
AIC_dd_1<-2*dd_1$df-2*dd_1$loglik
AIC_dd_2<-2*dd_2$df-2*dd_2$loglik
AIC_dd_2_1<-2*dd_2_1$df-2*dd_2_1$loglik
AIC_dd_2_2<-2*dd_2_2$df-2*dd_2_2$loglik
AIC_dd_3<-2*dd_3$df-2*dd_3$loglik
AIC_dd_4<-2*dd_4$df-2*dd_4$loglik
AIC_dd_4.1<-2*dd_4.1$df-2*dd_4.1$loglik
#AIC_dd_4.2<-2*dd_4.2$df-2*dd_4.2$loglik
AIC_dd_5<-2*dd_5$df-2*dd_5$loglik

AIC_dd_SR_1<-2*dd_SR_1$df-2*dd_SR_1$loglik
AIC_dd_SR_2<-2*dd_SR_2$df-2*dd_SR_2$loglik
AIC_dd_SR_2_1<-2*dd_SR_2_1$df-2*dd_SR_2_1$loglik
AIC_dd_SR_2_2<-2*dd_SR_2_2$df-2*dd_SR_2_2$loglik
AIC_dd_SR_3<-2*dd_SR_3$df-2*dd_SR_3$loglik
AIC_dd_SR_4<-2*dd_SR_4$df-2*dd_SR_4$loglik
AIC_dd_SR_4.1<-2*dd_SR_4.1$df-2*dd_SR_4.1$loglik
#AIC_dd_SR_4.2<-2*dd_SR_4.2$df-2*dd_SR_4.2$loglik

AICs<-c(AIC_bd_0,AIC_bd_1,AIC_bd_2,AIC_bd_3,AIC_dd_1,AIC_dd_2,AIC_dd_2_1,AIC_dd_2_2,AIC_dd_3,
        AIC_dd_4,AIC_dd_4.1,AIC_dd_5,AIC_dd_SR_1,AIC_dd_SR_2,AIC_dd_SR_2_1,AIC_dd_SR_2_2,AIC_dd_SR_3,AIC_dd_SR_4,AIC_dd_SR_4.1)

akaike.weights(AICs)
$deltaAIC
[1] 45.19130 45.19130 49.66640 49.62848 43.28770 42.65518 42.67986 37.94021 42.71269  0.00000 28.77553
[12] 44.71269 44.77055 45.22676 43.39211 61.53170 45.63449 47.79961 38.92981

$rel.LL
[1] 1.537563e-10 1.537563e-10 1.640888e-11 1.672291e-11 3.982856e-10 5.464431e-10 5.397425e-10
[8] 5.772827e-09 5.309552e-10 1.000000e+00 5.642528e-07 1.953273e-10 1.897570e-10 1.510543e-10
[15] 3.780275e-10 4.350712e-14 1.231961e-10 4.172978e-11 3.519651e-09

$weights
[1] 1.537563e-10 1.537563e-10 1.640887e-11 1.672290e-11 3.982854e-10 5.464428e-10 5.397422e-10
[8] 5.772823e-09 5.309549e-10 9.999994e-01 5.642525e-07 1.953271e-10 1.897569e-10 1.510542e-10
[15] 3.780273e-10 4.350709e-14 1.231960e-10 4.172975e-11 3.519649e-09

