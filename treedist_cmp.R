library(ggplot2)
spedist<-read.table("spetredist3.txt", header = T)
head(spedist)
ggplot(spedist, aes(x=reorder(name, TT), y=TT)) + geom_point(size = 3)  +  
  theme(axis.text.x = element_text(angle=60, hjust=1), panel.grid.major.y = element_blank()) +
  ylab ("Normalized Triple Distance") + xlab("tree pairs") + ylim(0, 0.08)

ggsave("spetredist3.tiff", dpi = 300)


dist<-read.csv("dist_sum.csv", header = T)
 
head(dist)
 
ggplot(dist, aes(x=median, fill=type)) +
  geom_histogram(position="identity", alpha=0.4, binwidth = 0.005)+scale_fill_discrete(name="gene-tree type",
                     breaks=c("bsdist", "ppdist"),
                     labels=c("BS tree", "PP tree")) +
          geom_vline(xintercept = 0.55, color = "coral3", linetype = "longdash") +
  
  geom_vline(xintercept = 0.42, color = "dodgerblue", linetype = "longdash") + 
  xlab("Normalized tree distances") + ylab("Density")

ggsave("bs_pp_dist.tiff", dpi = 300)
