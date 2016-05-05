library(ape)

###loading misfits script
source("misfitssingle.R")

###loading tree file
tree<-read.tree("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal2.txt")   

out<-misfitssingle (tree, filename="divmodels_misfits_cal2.txt")  ### you could replace the filename with whatever you what

out
