library(TreePar)
library(ape)
library(DDD)

tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)

dd_4_c2_r<-dd_ML(x,ddmodel=4,cond=2,res=5000)

save.image(file= "ppmpDDD_06.RData")
quit()