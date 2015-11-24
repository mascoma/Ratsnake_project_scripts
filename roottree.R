library(ape)
library(stringi)
files <- list.files(path = "/Users/Xin/Documents/Rat/mpest_1-2.5/mac/two_filters_pp/g6_BItrees", pattern="*.txt", full.names=T, recursive=F)
length(files)
stri_length(files[1])
stri_sub(files[1], 73, (stri_length(files[1])-4))
trees <- lapply(files, function(x) read.tree(x))
 
  for (i in 1 : length (trees)){
    tree1 <- root(trees[[i]], "I0111", resolve.root = T)
    genename <- stri_sub(files[i], 73, (stri_length(files[i])-4))
    name <- paste("g6BIroot_", genename, ".txt", sep = "")
    write.tree(tree1, file= name, append = T)
  }
 
files <- list.files(path = "/Users/Xin/Documents/Rat/mpest_1-2.5/mac/two_filters_pp/g5_pptreelines/g5_output", pattern="*.tre", full.names=T, recursive=F)
length(files)
for (i in 1:length(files)){
  tree <- read.nexus(files[i])
  write.tree(tree, file = "g5pp_mpest.txt", append = T)
}

trees<-read.tree("BI_trees_filtered.txt")
for (i in 1: length(trees)){
  tree1 <- root(trees[[i]], "I0111_FTB2500_Nerodia", resolve.root = T)
   
  write.tree(tree1, file= "BI_trees_filtered.txt", append = T) 
}