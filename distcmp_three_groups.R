library(ggplot2)

 
BItree_nonclk <-read.table("/Users/Xin/Documents/Rat/data_analysis/BItre_complete/BI_nonclkdist.txt", header = T)  
 
BItree_ppdist <-read.table("/Users/Xin/Documents/Rat/data_analysis/BItre_complete/BI_ppdist.txt", header = T)  

BItree_rest <-read.table("/Users/Xin/Documents/Rat/data_analysis/BItre_complete/BI_restdist.txt", header = T)  
head(BItree_nonclk)
head(BItree_ppdist)
head(BItree_rest)



  


nonclkqt<-BItree_nonclk[[5]]
quantile(nonclkqt)
mean(nonclkqt)
sd(nonclkqt)
median(nonclkqt)
range(nonclkqt)
plot(density(nonclkqt))
 
 

ppdistqt<-BItree_ppdist[[5]]
quantile(ppdistqt)
mean(ppdistqt)
sd(ppdistqt)
median(ppdistqt)
range(ppdistqt)
plot(density(ppdistqt))

 

restqt<-BItree_rest[[5]]
quantile(restqt)
mean(restqt)
sd(restqt)
median(restqt)
range(restqt)
plot(density(restqt))

 
BIdists <- matrix(data =NA, nrow = (length(BItree_nonclk[[5]]) + length(BItree_ppdist[[5]]) + length(BItree_rest[[5]])), ncol = 2)
BIdists[1:length(BItree_nonclk[[5]]), 1] <- 0
BIdists[(length(BItree_nonclk[[5]])+1):(length(BItree_nonclk[[5]]) + length(BItree_ppdist[[5]])), 1] <-1
BIdists[(length(BItree_nonclk[[5]])+length(BItree_ppdist[[5]])+1):(length(BItree_nonclk[[5]])+length(BItree_ppdist[[5]]) + length(BItree_rest[[5]])), 1] <-2
BIdists[1:length(BItree_nonclk[[5]]), 2] <-  nonclkqt
BIdists[(length(BItree_nonclk[[5]])+1):(length(BItree_nonclk[[5]]) + length(BItree_ppdist[[5]])), 2] <-ppdistqt
BIdists[(length(BItree_nonclk[[5]])+length(BItree_ppdist[[5]])+1):(length(BItree_nonclk[[5]])+length(BItree_ppdist[[5]]) + length(BItree_rest[[5]])), 2] <-restqt

BIdists<-as.data.frame(BIdists)
head(BIdists)
tail(BIdists)
colorder <- c( "khaki3", "darkseagreen", "lightpink3")
ggplot(BIdists, aes(x = factor(V1), y = V2)) +  geom_violin() +
geom_boxplot(width= 0.1, aes(fill=factor(V1))) + 
  scale_color_manual(breaks=colorder, # color scale (for points)
                     values=colorder,
                     labels=c("non-clock","large-variation", "rest"),
                     name="Group") +
  scale_fill_manual(breaks=colorder,  # fill scale (for boxes)
                    values=colorder,
                    labels=c("non-clock","large_variation", "rest"),
                    name="Group") +
  scale_x_discrete( labels=c("non-clock","large_variation", "rest")) + 
  ylab("Normalized Triple Distances") + xlab("")
ggsave("distcmp.tiff", dpi = 300)

 t.test(nonclkqt, restqt)

data:  nonclkqt and restqt
t = 8.8234, df = 3247.181, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.01860498 0.02923603
sample estimates:
  mean of x mean of y 
0.4395266 0.4156061 

t.test(ppdistqt, restqt)

data:  ppdistqt and restqt
t = 36.1832, df = 677.664, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.2415004 0.2692142
sample estimates:
  mean of x mean of y 
0.6709634 0.4156061 
 
#fisher's exact test
m1 = matrix(,2,2)
m1[1,] = c(8, 68)
m1[2,] = c(28, 199)
