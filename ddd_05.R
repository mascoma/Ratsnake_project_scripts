library(TreePar)
library(ape)
library(DDD)

tree<-read.tree("ppmpchrpltree.txt")
x<-getx(tree)

dd_4_c3_r<-dd_ML(x,ddmodel=4,res=5000)

save.image(file= "ppmpDDD_05.RData")
quit()