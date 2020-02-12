#library(plyr)
library(MSnbase)
library(reshape)
library(lattice)
library(ggplot2)
library(limma)
library(RColorBrewer)
library(ROTS)

library(seqinr)

library(ggrepel)
library(MKmisc) # glog2

library(gplots)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(gtools)

X11.options(width=5, height=5, xpos=1200, ypos=500)
options(editor="/usr/bin/vim")
options(stringsAsFactors=FALSE)


# LOAD DATA
#CHANGE input file name

Meta.df <- read.table("PCA_BTZ_NOT_comparison_Matched_noOutliers.csv", header=TRUE, row.names="MWRT", quote='\"', sep=',', comment.char='')
print(Meta.df[1:2])
meta.pca <- prcomp(Meta.df[,c(1:2)], center = TRUE, scale. = TRUE)

summary(meta.pca)
str(meta.pca)

p <- ggbiplot(meta.pca, labels=rownames(Meta.df)) +
  ggtitle('PCA of BTZ, NOT') +
  theme_classic()

tiff("PCA_BTZ_NOT_comparison_Matched_Labelled_PCA_noOutliers.tiff", width = 6, height = 4, units = 'in', res = 600)

plot(p)
dev.off()
