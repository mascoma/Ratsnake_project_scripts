library(ape)
library(ggtree)
bstree = read.nexus("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/bsmptree2.tre")
pptree = read.nexus("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/ppmptree2.tre")


#tips = c("Coelognathus_helena", "Coelognathis_subradiatus", "Coelognathus_erythrurus", "Coelognathus_flavolineatus", "Coelognathus_radiatus",
#         "Gonyosoma_oxycephalum", "Nerodia_sipedon", "Heterodon_platirhinos", "Gonyosoma_frenatum", "Gonyosoma_prasinum", "Gonyosoma_boulengeri", 
#         "Coluber_constrictor", "Drymobius_margaritiferus", "Gyalopion_canum", "Hapsidophrys_lineatus", "Hemorrhois_ravergieri","Tantilla_coronata")

#bs = drop.tip(bstree, tips)
#pp = drop.tip(pptree, tips)

#bs =  chronopl(bstree, lambda = 0.01, age.min = 10)
#pp = chronopl(pptree, lambda = 0.01, age.min = 10)
par(mar = c(0, 0, 0, 0))
og = c(1:4, 144:172)
colo = rep("black", Nedge(pptree))
colo[og] = "grey50"
png("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/ppmptree2.png", width = 1300, height = 1872)
plot(ladderize(pptree), adj = 0.5, cex = 1.5, x.lim = 30, edge.width = 3, edge.color = colo,font = 4, label.offset = 0.2)
pptree$node.label = c("-", 1, 0.54, 0.93, 0.97, 0.98, 0.98, 0.8, 0.98, 0.98, 0.99, 0.9, 0.78, 0.58, 0.69, 0.81, 0.52, 0.96,0.41, 0.49, 
                      0.51, 0.43, 0.58, 0.58, 0.7, 0.95, 0.97, 0.55, 0.86, 0.64, 0.94, 0.99, 0.23, 0.99, 0.45, 0.88, 0.62, 0.74, 0.61, 0.96, 
                      0.9, 0.53, 0.92, 0.93, 0.96, 0.81, 0.83, 0.86, 0.88, 0.96, 0.97, 0.99, 1, 0.96, 1, 0.8, 0.99, 1, 1, 0.78, 0.96, 1, 
                      0.9, 0.92, 0.92, 0.9, 0.92, 0.53, 0.99, 0.97, 0.56, 0.98, 1, 0.92, 0.74, 0.98, 0.91,0.98, 0.99,  0.95, 0.98, 0.99, 0.93, 0.98,1, 0.99)
#nodelabels()
#edgelabels()
nodecolor = rep("black", length(pptree$node.label))
nodecolor[1:17] = "grey50"
nodelabels(pptree$node.label, frame = "none", cex = 1.4, adj = c(1.05, -0.5), col = nodecolor, font = 2)
dev.off()

par(mar = c(0, 0, 0, 0))
og = c(1:3, 143:172)
colo = rep("black", Nedge(bstree))
colo[og] = "grey50"
png("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/bsmptree2.png", width = 1070, height = 1872)
plot(ladderize(bstree),  show.tip.label = F,  direction = "l" , edge.width = 3, edge.color = colo)
#nodelabels()
bstree$node.label = c("-", 1,0.58, 0.84, 0.92, 0.94, 0.96, 0.36, 0.93, 0.88, 0.96, 0.72, 0.63, 0.54, 0.67, 0.75, 0.8, 0.87, 0.54, 0.61, 0.46, 
                      0.48, 0.61, 0.49, 0.54, 0.88, 0.92, 0.41, 0.67, 0.56, 0.71, 0.95, 0.18, 0.96, 0.38, 0.79, 0.39, 0.37, 0.37, 0.91, 0.68, 0.46,
                      0.74, 0.77, 0.9, 0.8, 0.9, 0.33, 0.38, 0.76, 0.64, 0.9, 0.97, 0.73, 0.98, 0.62, 0.75, 0.94, 0.96, 0.74, 0.67, 0.97, 0.5,
                      0.46, 0.54, 0.4, 0.62, 0.44, 0.97, 0.9, 0.74, 0.83, 0.98, 0.72, 0.54, 0.98, 0.56, 0.91, 0.98, 0.75, 0.93, 0.98, 0.7, 0.91, 0.97, 0.96)
#edgelabels()
nodecolor = rep("black", length(bstree$node.label))
nodecolor[1:17] = "grey50"
nodelabels(bstree$node.label, frame = "none", cex = 1.4, adj = c(-0.09, -0.5), col = nodecolor, font = 2)
dev.off()