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

TMTData.df <- read.table("200113L_SAM6526_BY53_CC_PO4_pH_All_KB_Percolator_sitelevel.csv", header=TRUE,quote='\"', sep=',', comment.char='')
tmt <- as.matrix(TMTData.df[16:23])

TMTImputed.df  <- read.table("200113L_SAM6526_BY53_CC_PO4_pH_All_KB_Percolator_sitelevel_Imputed.csv", header=TRUE,quote='\"', sep=',', comment.char='')
tmtimp <- as.matrix(TMTImputed.df[16:23])

y <- normalizeVSN(tmt)
ydata <- melt(tmt) 

x <- normalizeVSN(tmtimp)
xdata <- melt(tmtimp)

p <- ggplot(ydata, aes(X2, log10(value))) +
  geom_boxplot() + stat_boxplot(geom = 'errorbar') +
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1))
tiff("vsnNormalizationBox_notImputed.tiff", width = 4, height = 4, units = 'in', res=600)
plot(p)
dev.off()

blankttest.df  <- read.table("Imputed_vs_NotImputed_vs_MCARimputed.csv", header=TRUE,quote='\"', sep=',', comment.char='')

p <- ggplot(blankttest.df, aes(x=LOG2FC, y=LOG10TTEST)) + 
  geom_point(aes(colour=GROUP),size = 1.5) +
  geom_hline(yintercept=-log(0.05,10)) +
  theme_classic()
p <- p + labs(title='Morphine vs Control Comparison of Imputed vs Not Imputed vs MCAR Imputed', x='log2(FC)', y='-log10(Nominal P-value)')
#change file name
tiff("Missing vs Blank vs Half_Volcano.tiff", width = 8, height = 6, units = 'in', res = 600)

plot(p)
dev.off()


TMTData.df <- read.table("PCA_Imputed.csv", header=TRUE,quote='\"', sep=',', comment.char='')
tmt <- as.matrix(TMTData.df)
#change file name
pd <- c('MOR', 'CTL')
pca <- prcomp(tmt, scale=FALSE)
df <- as.data.frame(pca$rotation[, 1:4])
df <- namerows(df, col.name='Samples')

df$Samples <- c('MOR', 'MOR', 'MOR', 'MOR', 'CTL', 'CTL', 'CTL', 'CTL')

p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2)
tiff("PCA_Imputed.tiff", width = 6, height = 4, units = 'in', res=600)
plot(p)
dev.off()
