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
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(gtools)

library(rgl)
library(pca3d)

X11.options(width=5, height=5, xpos=1200, ypos=500)
options(editor="/usr/bin/vim")
options(stringsAsFactors=FALSE)


# LOAD DATA
#CHANGE input file name

Meta.df <- read.table("VSN_normalized_noOutliers_format.csv", header=TRUE,quote='\"', sep=',', comment.char='')
print(Meta.df[1:3])
mx <- as.matrix(Meta.df[2:27])

pca <- prcomp(mx, scale=FALSE)
df <- as.data.frame(pca$rotation[, 1:4])
df <- namerows(df, col.name='Samples')
write.table(df,file="editing.csv")
MetaEDIT.df <- read.table("editing.csv", header=TRUE,quote='\"', sep=' ', comment.char='')

p <- ggplot(MetaEDIT.df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2)
tiff("PCA_nofileSep.tiff", width = 6, height = 4, units = 'in', res=600)
plot(p)
dev.off()
