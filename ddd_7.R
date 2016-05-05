library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test

dd_SR_2_2<-dd_SR_ML(x,ddmodel=2.2)
dd_SR_3<-dd_SR_ML(x,ddmodel=3)

save.image(file= "ppmpDDD7.RData")
quit()