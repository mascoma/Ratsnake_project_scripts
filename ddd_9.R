library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test


dd_SR_4.2<-dd_SR_ML(x,ddmodel=4.2)
save.image(file= "ppmpDDD9.RData")
quit()