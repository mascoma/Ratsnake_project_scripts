library(TreePar)
library(ape)
library(DDD)
tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)
dd_4_c1_r<-dd_ML(x,ddmodel=4,cond=0,res=5000)

save.image(file= "ppmpDDD_08.RData")
quit()