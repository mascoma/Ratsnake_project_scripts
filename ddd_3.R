library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test
dd_3<-dd_ML(x,ddmodel=3)
dd_4<-dd_ML(x,ddmodel=4)
save.image(file= "ppmpDDD3.RData")
quit()