library(TreePar)
library(ape)
library(DDD)

tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)

dd_4_c2<-dd_ML(x,ddmodel=4,cond=2)

save.image(file= "ppmpDDD_04.RData")
quit()