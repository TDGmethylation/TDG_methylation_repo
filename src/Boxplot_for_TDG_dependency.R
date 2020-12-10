library(gplots)
library(ggplot2)
library(ggpubr)

TDG_context_CpN <- read.csv("/Users/gerogehou/Desktop/R/TDG_context_analysis/TDG_FI_NNT_ACGT_NNN.tsv", 
                            sep = "\t", header = FALSE)

TDG_context_TAGN <- read.csv("/Users/gerogehou/Desktop/R/TDG_context_analysis/TDG_FI_NNTAG_ACGT_N.tsv", 
                             sep = "\t", header = FALSE)
TDG_context_TGCN <- read.csv("/Users/gerogehou/Desktop/R/TDG_context_analysis/TDG_FI_NNTGC_ACGT_N.tsv",
                             sep = "\t", header = FALSE)

TN <- ggplot(TDG_context_CpN, aes(x = V1, y = V2)) + geom_boxplot()
TAG <- ggplot(TDG_context_TAGN, aes(x = V1, y = V2)) + geom_boxplot()
TGC <- ggplot(TDG_context_TGCN, aes(x = V1, y = V2)) + geom_boxplot()

summary(TDG_context_CpN)
summary(TDG_context_TAGN$V2)
summary(TDG_context_TGCN$V2)

#mean(TDG_context_TGCN$V2[TDG_context_TGCN$V1[ ] == "A"])

##########################HEAT_MAP_FOR_FIGURE_3C######################################
R_methylation <- read.csv("/Volumes/Genome_data/Methylation_R_summary/R_for_all_cellline.tsv", 
                            sep = "\t", header = TRUE)
rownames(R_methylation) <- c("CpA", "CpC", "CpG", "CpT")

N_breaks = seq(-1, 1, length.out = 1000)
grad_color1 = colorpanel( sum( N_breaks[-1]<= 0 ), "blue", "white")
grad_color2 = colorpanel( sum( N_breaks[-1]> 0 ), "white", "red")
hm.colors = c(grad_color1, grad_color2)

heatmap.2(as.matrix(R_methylation),scale="none",breaks=N_breaks,col=hm.colors, Colv=FALSE,
          dendrogram="row",trace="none")