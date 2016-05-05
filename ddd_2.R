library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test
dd_2_1<-dd_ML(x,ddmodel=2.1)
dd_2_2<-dd_ML(x,ddmodel=2.2)
save.image(file= "ppmpDDD2.RData")
quit()
