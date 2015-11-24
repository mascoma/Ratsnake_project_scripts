 
library(ape)
library(stringi)
 
 
  trees <- read.tree("BItrees_complete.txt")
  for (j in 1 : length(trees)){
    tree1 <- root(trees[[j]], "I0111", resolve.root = T)      
    write.tree(tree1, file= "BIrooted_304.txt", append = T)
  }  
 


files <- list.files(path = "/Users/Xin/Documents/Rat/304geneloci/mb_analysis/subsampled/g6_BItrees", pattern="*.txt", full.names=T, recursive=F)
for (i in 1 : length(files)) {
  tree <- read.tree(files[i])
  treename <- paste("g6BIroot_gene", substr(files[i], 81, stri_length(files[i])),sep="")
 
    tree1 <- root(tree, "I0111", resolve.root = T)      
    write.tree(tree1, file= treename)
}

 

dir()
files <- list.files(path = "/Users/Xin/Desktop/mpest_1-2.5/mac/g5_pptrees", pattern="*.txt", full.names=T, recursive=F)
length(files)


splitbsgenetree<-function(files,n){
  ## files: listed all the tree file names; n: number of the trees in each tree file
  
  trees<-lapply(files, function(x) read.nexus(x)) 
  
  
  name<-vector("character",length=n)
  name2<-vector("character",length=n)
  for (a in 1 : n){
    
    name[a]<-paste("g5_pptreeline",a,sep="_")
    name2[a] <- paste(name[a], ".txt", sep="")
  }
  
  
  
  for (i in 1 : length(files)){
    
    tmp<-trees[[i]]
    for (j in 1 : n){
      write.tree(tmp[[j]],file=name2[j],append=T)
      
    }
    
  }
  
}

splitbsgenetree (files,1)
