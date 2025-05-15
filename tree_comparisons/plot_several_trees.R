### Plots several tree comparisons

library(ggplot2)
library(phytools)

setwd("~/tree_comparisons/")

root <- as.multiPhylo(read.tree("all.phyluce.trees.tre"))
# tree order: ast.68p.bs10, ast.68p, ast.total.bs10, ast.total, concat.68p, concat.total

# plot all trees overlapping (from Fig. 2)
pdf("incongruence_phyluce_astral_rooted_twoOptions.pdf", width = 10, height = 10)
densityTree(root,use.edge.length=FALSE,type="phylogram",nodes="inner", color = "purple", fsize = 0.7)
densityTree(root,use.edge.length=FALSE,type="cladogram",nodes="inner", color = "purple", fsize = 0.7)
dev.off()

### tanglegrams (supplmentals)

obj <- cophylo(root[[1]], root[[2]])

pdf("tangle_astral68.pdf", height = 10, width = 14)
plot(obj,type=c("phylogram","phylogram"),fsize=0.8,part=0.48,
     mar=c(0.1,0.1,2.1,0.1),pts=FALSE,lwd=2,link.type="curved")
left<-get("last_plot.cophylo",envir=.PlotPhyloEnv)[[1]]
ind<-2:Nnode(obj$trees[[1]])+Ntip(obj$trees[[1]])
points(left$xx[ind],left$yy[ind],
       cex=2*as.numeric(obj$trees[[1]]$node.label[2:Nnode(obj$trees[[1]])]),
       pch=16,col=make.transparent(palette()[3],0.75))
right<-get("last_plot.cophylo",envir=.PlotPhyloEnv)[[2]]
ind<-2:Nnode(obj$trees[[2]])+Ntip(obj$trees[[2]])
points(right$xx[ind],right$yy[ind],
       cex=2*as.numeric(obj$trees[[2]]$node.label[2:Nnode(obj$trees[[2]])]),
       pch=16,col=make.transparent(palette()[3],0.75))
legend("topleft",
       legend=paste(seq(20,100,by=20),"%",sep=""),
       pt.cex=seq(0.4,2,by=0.4),
       pch=16,col=make.transparent(palette()[3],0.75),
       title="bootstrap support",bty="n")
mtext("ASTRAL 68p",at=-0.5,adj=0)
mtext("ASTRAL 68p clades collapsed",at=0.02,adj=0)
dev.off()


### more tangles

obj1 <- cophylo(root[[1]], root[[2]])
obj2 <- cophylo(root[[3]], root[[4]])
obj3 <- cophylo(root[[5]], root[[6]])
obj4 <- cophylo(root[[1]], root[[6]])

pdf("tangle1.pdf", width = 10, height = 14)
plot(obj1,link.type="curved",link.lwd=3,link.lty="solid", link.col=make.transparent("purple",0.25),fsize=1)
dev.off()

pdf("tangle2.pdf", width = 10, height = 14)
plot(obj2,link.type="curved",link.lwd=3,link.lty="solid", link.col=make.transparent("purple",0.25),fsize=1)
dev.off()

pdf("tangle3.pdf", width = 10, height = 14)
plot(obj3,link.type="curved",link.lwd=3,link.lty="solid", link.col=make.transparent("purple",0.25),fsize=1)
dev.off()

pdf("tangle4.pdf", width = 16, height = 14)
plot(obj4,link.type="curved",link.lwd=3,link.lty="solid", link.col=make.transparent("purple",0.25),fsize=1)
dev.off()



### MDS of RF distance (phyluce only, shown in Fig 2)
matr <- read.csv("treedist.txt", sep = "\t", header = TRUE, row.names = 1)


makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

distance <- makeSymm(matr)

cols <- c("blue", "blue", "blue", "blue", "red", "red")
shapes <- c(0, 15,1, 16, 15, 16)

pdf("mds_trees.pdf", height = 8, width = 8)
plot(cmdscale(distance), xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", col = cols, pch = as.numeric(shapes))
text(cmdscale(distance),  cex=1, pos=3)
dev.off()

### plotting support values (not used)

##extract support values
##order: 68p, 68p bs10, total, total bs10
s1 <- root[[1]]$node.label
s2 <- root[[2]]$node.label
s3 <- root[[3]]$node.label
s4 <- root[[4]]$node.label
s5 <- new[[1]]$node.label
s6 <- new[[2]]$node.label

df <- data.frame(s1, s2, s3, s4, s5, s6)
colnames(df) <- c("68p", "68p bs10", "total", "total bs10", "concat 68p", "concat total")
write.csv(df, "suppvalues.csv")

supp <- read.csv("suppvalues.csv")

pdf("supp.violin.pdf", height = 7, width = 10)
ggplot(supp, aes(x = treatment, y = support.values, fill = treatment)) + geom_violin() + 
  scale_fill_brewer(palette = "Set2", alpha(0.70)) + theme(legend.position="none") + labs(x="Type of analysis", y = "Support values") + 
  geom_boxplot(width=0.1, color="black", alpha=0.6, lwd = 0.2)
dev.off()


