library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ingrouptree.txt")
x<-getx(tree)
###DDD test

dd_4_c0<-dd_ML(x,ddmodel=4,cond=0)
save.image(file= "ppmpDDD_02.RData")
quit()

