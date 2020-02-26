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

library(rgl)
library(pca3d)
X11.options(width=5, height=5, xpos=1200, ypos=500)
options(editor="/usr/bin/vim")
options(stringsAsFactors=FALSE)


# LOAD DATA
#CHANGE input file name

Meta.df <- read.table("CellMatrix_LEN_LBH_BTZ_NOTFilledChannels.csv", header=TRUE,quote='\"', sep=',', comment.char='')
mx <- as.matrix(Meta.df[2:38])
meta.pca <- prcomp(Meta.df[,c(1:3)], center = TRUE, scale. = TRUE)

mdata <- melt(Meta.df, id=c("Row"))


y <- normalizeVSN(mx)
ydata <- melt(y)
write.table(y,file="VSN_normalized.csv")

xNoOutl <- read.table("VSN_normalized_noOutliers.csv", header=TRUE,quote='\"', sep=',', comment.char='')
xdata <- melt(xNoOutl)
p <- ggplot(xdata, aes(variable, log10(value))) +
  geom_boxplot() + stat_boxplot(geom = 'errorbar') +
theme(axis.text.x = element_text(angle=90, hjust=1))
tiff("vsnNormalizationBox_noOutliers.tiff", width = 4, height = 4, units = 'in', res=600)
plot(p)
dev.off()
