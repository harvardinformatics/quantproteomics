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

TMTFemale.df  <- read.table("200113L_SAM6526_BY53_CC_PO4_pH_All_KB_Percolator_sitelevel_Female.csv", header=TRUE,quote='\"', sep=',', comment.char='')

TMTMale.df  <- read.table("200113L_SAM6526_BY53_CC_PO4_pH_All_KB_Percolator_sitelevel_Male.csv", header=TRUE,quote='\"', sep=',', comment.char='')

y <- normalizeVSN(tmt)
ydata <- melt(tmt) 

x <- normalizeVSN(tmtimp)
xdata <- melt(tmtimp)

p <- ggplot(xdata, aes(X2, log10(value))) +
  geom_boxplot() + stat_boxplot(geom = 'errorbar') +
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1))
tiff("vsnNormalizationBox_Imputed.tiff", width = 4, height = 4, units = 'in', res=200)
plot(p)
dev.off()

blankttest.df  <- read.table("Imputed_vs_NotImputed_vs_MCARimputed.csv", header=TRUE,quote='\"', sep=',', comment.char='')
blanktNotImputed.index <- grep("^NotImputed", blankttest.df$GROUP)
blanktNotImputed.df <- blankttest.df[blanktNotImputed.index,]

blanktImputed.index <- grep("^Imputed", blankttest.df$GROUP)
blanktImputed.df <- blankttest.df[blanktImputed.index,]

blanktMCARImputed.index <- grep("^MCARImputed", blankttest.df$GROUP)
blanktMCARImputed.df <- blankttest.df[blanktMCARImputed.index,]

twoBytwo.df  <- read.table("MCAR_NotImputed_2x2.csv", header=TRUE,quote='\"', sep=',', comment.char='')
notImputed.index <- grep("NotImputed", twoBytwo.df$GROUP)
notImputed.df <- twoBytwo.df[notImputed.index,]

MCARImputed.index <- grep("MCAR", twoBytwo.df$GROUP)
MCARImputed.df <- twoBytwo.df[MCARImputed.index,]

p <- ggplot(notImputed.df, aes(x=LOG2FC, y=LOG10TTEST)) + 
  geom_point(color='red', size = 1.5) +
  geom_hline(yintercept=-log(0.05,10)) +
  theme_classic()
p <- p + labs(title='Morphine vs Control Not Imputed', x='log2(FC)', y='-log10(Nominal P-value)')
#change file name
tiff("NOTImputed_Volcano.tiff", width = 8, height = 6, units = 'in', res = 200)

plot(p)
dev.off()


p <- ggplot(blanktNotImputed.df, aes(x=BLANKS, y=LOG2FC)) + 
  geom_point(color='red',size = 1.5) +
  theme_classic()
p <- p + labs(title='Morphine vs Control Comparison of logFC vs Number of Blanks', x='Number of Blanks', y='log2FC')
#change file name
tiff("FCBlanks_NotImputed.tiff", width = 8, height = 6, units = 'in', res = 200)

plot(p)
dev.off()

p <- ggplot(TMTFemale.df, aes(x=R1, y=R2)) + 
  geom_point(color='red',size = 1.5) +
  theme_classic()
p <- p + labs(title='Morphine vs Control Ratio 1 and Ratio 2 Female', x='T1/C1', y='T2/C2')
#change file name
tiff("R1_R2_pointFemale.tiff", width = 8, height = 6, units = 'in', res = 200)

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
