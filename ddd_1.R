library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
###DDD test
bd_0<-bd_ML(x,tdmodel=0) 
bd_1<-bd_ML(x,tdmodel=1)
bd_2<-bd_ML(x,tdmodel=2)
bd_3<-bd_ML(x,tdmodel=3)

dd_1<-dd_ML(x,ddmodel=1)
dd_2<-dd_ML(x,ddmodel=2)
save.image(file= "ppmpDDD1.RData")
quit()