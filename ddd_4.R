library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test
dd_4.1<-dd_ML(x,ddmodel=4.1)
dd_4.2<-dd_ML(x,ddmodel=4.2)
save.image(file= "ppmpDDD4.RData")
quit()