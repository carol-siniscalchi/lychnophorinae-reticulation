library(ggmagnify)
library(ggpubr)
library(phytools)
library(ggplot2)
library(plyr)

### compiling dataset

data <- read.csv("~/assembly_comparisons/hybpiper.plots.data.tsv", sep = "\t")
phyd <- read.csv("~/assembly_comparisons/phyluce.total.stats.csv")

total <- subset(data, treatment == "withParalogs")
para <- subset(data, treatment == "noParalogs")



comparison <- rbind(phyd, para)

summary(total)
summary(para)

ggplot(total, aes(x=Alignment_length)) +
        geom_histogram()

ggplot(para, aes(x=Alignment_length)) +
  geom_histogram()

mu <- ddply(data, "treatment", summarise, grp.mean=mean(Alignment_length))

ggplot(data, aes(x=Alignment_length, fill = treatment, color = treatment)) +
  geom_histogram(alpha=0.5)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=treatment),
            linetype="dashed")


ggplot(phyd, aes(x=Alignment_length)) +
  geom_histogram()

mu2 <- ddply(comparison, "treatment", summarise, grp.mean=mean(Alignment_length))

### comparing phyluce and hybpiper without paralogs
ggplot(comparison, aes(x=Alignment_length, fill = treatment, color = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5)+
  geom_vline(data=mu2, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed")+
  xlim(0, 1000)

mu3 <- ddply(comparison, "treatment", summarise, grp.mean=mean(Missing_percent))

## Missing %
mis1 <- ggplot(comparison, aes(x=Missing_percent, fill = treatment)) +
  geom_histogram(binwidth = 5, alpha=0.5, color = "gray30")+
  labs(y="Number of loci", x = "Percentage of Missing Data") +
  scale_fill_discrete(name = "Assembly", labels = c("HybPiper", "phyluce")) +
  geom_vline(data=mu3, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed", show.legend = F)

mis2 <- ggplot(comparison, aes(y=Missing_percent, fill = treatment)) +
  geom_boxplot(color = "gray30", show.legend = F) +
  labs(y="Percentage of missing data", x = "")


pdf("comparison_missing.pdf", width = 12, height = 6)
ggarrange(mis1, mis2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()

## alignment length
ali1 <- ggplot(comparison, aes(x=Alignment_length, fill = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5, color = "gray30", show.legend = F)+
  geom_vline(data=mu2, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed", , show.legend = F) +
  labs(y="Number of loci", x = "Alignment length (bp)") + 
  scale_fill_discrete(name = "Assembly", labels = c("HybPiper", "phyluce")) 

ali2 <- ggplot(comparison, aes(y=Alignment_length, fill = treatment)) +
  geom_boxplot(color = "gray30", show.legend = F) +
  labs(y="Alignment length (bp)", x = "")

pdf("comparison_length.pdf", width = 15, height = 6)
ggarrange(ali1, ali2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()

pdf("inset.pdf", width = 8, height = 4)
ggplot(comparison, aes(x=Alignment_length, fill = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5, color = "gray30", show.legend = F)+
  geom_vline(data=mu2, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed", , show.legend = F) +
  labs(y="Number of loci", x = "Alignment length (bp)") + 
  scale_fill_discrete(name = "Assembly", labels = c("HybPiper", "phyluce")) + xlim(0, 1000)
dev.off()


### comparing phyluce total vs. 68p

phy68 <- read.csv("phyluce.68p.stats.csv")
colnames(phy68)[1] <- "treatment"

phycomp <- rbind(phyd, phy68)
mu4 <- ddply(phycomp, "treatment", summarise, grp.mean=mean(Alignment_length))

p1 <- ggplot(phycomp, aes(x=Alignment_length, fill = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5, color = "gray30")+
  geom_vline(data=mu4, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed", , show.legend = F) +
  labs(y="Number of loci", x = "Alignment length (bp)") + 
  scale_fill_manual(values = c("#00A6A6", "#BBDEF0")) +
  theme(legend.position="top")

p2 <- ggplot(phycomp, aes(y=Alignment_length, fill = treatment)) +
  geom_boxplot(color = "gray30") +
  scale_fill_manual(values = c("#00A6A6", "#BBDEF0")) +
  labs(y="Alignment length (bp)", x = "") +
  theme(legend.position="top")

pdf("comparison_phyluce.pdf", width = 12, height = 5)
ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()

### HybPiper

hybdata <- read.csv("hybpiper.noparalogs.75total.stats.csv")
mu5 <- ddply(comp, "treatment", summarise, grp.mean=mean(Alignment_length))

sums <- aggregate( hybdata , by=list(hybdata$treatment) , FUN=summary)
write.csv(sums, "summary_hp.csv")

partial <- subset(hybdata, treatment == "partial_noParalogs")
total <- subset(hybdata, treatment == "total_noParalogs")
comp <- rbind(para, partial)

p1 <- ggplot(comp, aes(x=Alignment_length, fill = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5, color = "gray30")+
  labs(y="Number of loci", x = "Alignment length (bp)") + 
  scale_fill_manual(values = c("#00A6A6", "#BBDEF0")) +
  theme(legend.position="top")+ 
  geom_vline(data=mu5, aes(xintercept=grp.mean, color=treatment), linetype="dashed", show.legend = F)

ggplot(comp, aes(x=Alignment_length, fill = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5, color = "gray30")+
  labs(y="Number of loci", x = "Alignment length (bp)") + 
  theme(legend.position="top")+ 
  geom_vline(data=mu5, aes(xintercept=grp.mean, color=treatment), linetype="dashed", show.legend = F)

p2 <- ggplot(comp, aes(y=Alignment_length, fill = treatment)) +
  geom_boxplot(color = "gray30") +
  scale_fill_manual(values = c("#00A6A6", "#BBDEF0")) +
  labs(y="Alignment length (bp)", x = "") +
  theme(legend.position="top")

pdf("comparison_hybpiper.pdf", width = 12, height = 5)
ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()

reduced <- rbind(phy68, partial)
mu6 <- ddply(reduced, "treatment", summarise, grp.mean=mean(Alignment_length))

r1 <- ggplot(reduced, aes(x=Alignment_length, fill = treatment)) +
  geom_histogram(binwidth = 50, alpha=0.5, color = "gray30")+
  labs(y="Number of loci", x = "Alignment length (bp)") + 
  scale_fill_manual(values = c("#00A6A6", "#BBDEF0")) +
  theme(legend.position="top")+ 
  geom_vline(data=mu6, aes(xintercept=grp.mean, color=treatment), linetype="dashed", show.legend = F)

r2 <- ggplot(reduced, aes(y=Alignment_length, fill = treatment)) +
  geom_boxplot(color = "gray30") +
  scale_fill_manual(values = c("#00A6A6", "#BBDEF0")) +
  labs(y="Alignment length (bp)", x = "") +
  theme(legend.position="top")


pdf("comparison_hybpiper75-phyluce68.pdf", width = 12, height = 5)
ggarrange(r1, r2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()

### RF and MDS

all_trees <- read.tree("~/assembly_comparisons/all.wHybPiper.tre")

multiRF(all_trees)

distance <- read.csv("~/assembly_comparisons/new_treedist.txt", sep = "\t", header = TRUE, row.names = 1)

cols <- c("blue", "blue", "blue", "blue", "red", "red", "cyan", "cyan")
shapes <- c(0, 15, 1, 16, 15, 16, 17, 17)

pdf("mds_trees_withHyb_true.pdf", height = 8, width = 8)
plot(cmdscale(distance), xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", col = cols, pch = as.numeric(shapes))
text(cmdscale(distance),  cex=1, pos=3)
dev.off()

#### paralogs between phyluce and HybPiper ###

parag <- read.csv("~/assembly_comparisons/paralogs_long.csv", header = TRUE)

pdf("paralogs_comparison.pdf", width = 8, height = 4)
ggplot(parag, aes(x=Treatment, y=Genes, fill=Classification)) + 
  geom_violin(position = position_dodge(0.9), alpha = 0.6) +
  scale_fill_manual(values = c("dodgerblue", "firebrick")) +
  geom_boxplot(width=0.05, position = position_dodge(0.9))
dev.off()


