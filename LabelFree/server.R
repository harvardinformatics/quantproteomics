options(shiny.maxRequestSize=300*1024^2) 

library(mzR)
library(mwshiny)
library(stringr)

shinyServer(function(input, output, session) {

  
  observeEvent(input$extractmzml, {

  mz <- openMSfile(input$entermzml)
  fileName(mz)
  print(instrumentInfo(mz))
  
  getScan <- function(msScan, scanNum) {
    getData <- peaks(msScan, scan=scanNum)
    
  }
  PSMfile <- data.frame(read.csv(input$enterpsm))
  scannums <- PSMfile$First.Scan
  mzmlfiles <- PSMfile$Spectrum.File
  for (i in 1:nrow(PSMfile)){
    sfile <- str_replace(PSMfile[i,]$Spectrum.File, '.raw', '.mzML')
    MSdata <- openMSfile(sfile)
    x<-getScan(MSdata, PSMfile[i,]$First.Scan)
    sfile_ext <- str_replace(PSMfile[i,]$Spectrum.File, '.raw', '')
    write.csv(x,paste('spectra/compound_spectra ',sfile_ext,PSMfile[i,]$First.Scan,'.csv'), row.names = F)
  }
  })
  
  observeEvent(input$calcperc, {
    
    system(paste("python3 CalcAbundancePerc.py", input$enterpsm, wait=FALSE))
    
  })
  
})