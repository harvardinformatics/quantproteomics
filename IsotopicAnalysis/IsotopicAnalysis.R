library(mzR)

fileMS <- system.file()
MSdata <- openMSfile('SUB09034_SPL007.mzML')
fileName(MSdata)
y <- runInfo(MSdata)

getScan <- function(msScan, scanNum, retentionTime) {
  getData<-peaks(msScan, scan=scanNum)
  
}

x<-getScan(MSdata, 4, 71.45420)
x
#x[1,][1]
write.table(y,'SUB09034_SPL007_info.txt', col.names = F, row.names = F)
write.table(x,'SUB09034_SPL007.txt', col.names = F, row.names = F)
