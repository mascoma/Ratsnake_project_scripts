library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test

dd_5<-dd_ML(x,ddmodel=5)

dd_SR_1<-dd_SR_ML(x,ddmodel=1)
save.image(file= "ppmpDDD5.RData")
quit()