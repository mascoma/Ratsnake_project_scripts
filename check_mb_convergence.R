pstat<-read.table("pstatall.txt", sep="\t", header=T)
head(pstat)
min(pstat[,"minESS"])
range(pstat[, "PSRF"])
