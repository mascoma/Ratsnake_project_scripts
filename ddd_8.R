library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test



dd_SR_4<-dd_SR_ML(x,ddmodel=4)
dd_SR_4.1<-dd_SR_ML(x,ddmodel=4.1)

save.image(file= "ppmpDDD8.RData")
quit()