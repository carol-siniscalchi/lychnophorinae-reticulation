library(ggplot2)

long <- read.csv("loglik_long.csv", header = TRUE)

colnames(long) <- c("H_number", "Ploglik", "Analysis")

pdf(file = "loglik.pdf", width = 10, height = 10)
ggplot(long, aes(x=H_number, y=Ploglik, group=Analysis, color=Analysis)) + 
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set1")
dev.off()