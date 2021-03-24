library(mzR)

fileMS <- system.file()
MSdata <- openMSfile('SUB09034_SPL007.mzML')
fileName(MSdata)
runInfo(MSdata)

head(peaks(MSdata))



#write.table(head(peaks(MSdata)),'SUB09034_SPL007.txt')
