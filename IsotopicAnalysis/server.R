options(shiny.maxRequestSize=300*1024^2) 

library(mzR)
library(mwshiny)
library(dplyr)
shinyServer(function(input, output, session) {
  observeEvent(input$calciso, {
    isofor <- as.character(input$isoformula)
    print(isofor)
    system(paste("python3 IsotopeAnalysis.py ", isofor, " ", input$ionmode, " ", input$adduct, wait=FALSE))
    
  })
  
  observeEvent(input$extractmzml, {
  fileMS <- system.file()
  MSdata <- openMSfile(input$entermzml)
  fileName(MSdata)
  y <- runInfo(MSdata)
  
  getScan <- function(msScan, scanNum) {
    getData <- peaks(msScan, scan=scanNum)
    
  }
  
  x<-getScan(MSdata, as.numeric(input$scannum))
  x
  write.csv(x,'compound_spectra.csv', row.names = F)
  
  })
  
  observeEvent(input$calcperc, {
    
    system(paste("python3 CalcAbundancePerc.py compound_spectra.csv mz_calc.csv", wait=FALSE))
    
  })
  
})