library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test



dd_SR_2<-dd_SR_ML(x,ddmodel=2)
dd_SR_2_1<-dd_SR_ML(x,ddmodel=2.1)
save.image(file= "ppmpDDD6.RData")
quit()