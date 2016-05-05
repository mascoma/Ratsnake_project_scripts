library(TreePar)
library(ape)
library(DDD)
tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)

dd_4_c1<-dd_ML(x,ddmodel=4,cond=1)
 
save.image(file= "ppmpDDD_03.RData")
quit()
