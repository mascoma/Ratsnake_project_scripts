library(TreePar)
library(ape)
library(DDD)


tree<-read.tree("ingrouptree.txt")
x<-getx(tree)
###DDD test
bd_0<-bd_ML(x,tdmodel=0) 
bd_1<-bd_ML(x,tdmodel=1)
bd_2<-bd_ML(x,tdmodel=2)
bd_3<-bd_ML(x,tdmodel=3)

bd_0_c0<-bd_ML(x,tdmodel=0, cond=0) 
bd_1_c0<-bd_ML(x,tdmodel=1,cond=0)
bd_2_c0<-bd_ML(x,tdmodel=2,cond=0)
bd_3_c0<-bd_ML(x,tdmodel=3,cond=0)



save.image(file= "ppmpDDD_01.RData")
quit()