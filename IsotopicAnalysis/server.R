options(shiny.maxRequestSize=300*1024^2) 

library(mzR)
library(mwshiny)

shinyServer(function(input, output, session) {
  observeEvent(input$calciso, {
    isofor <- as.character(input$isoformula)
    print(isofor)
    system(paste("python3 IsotopeAnalysis.py ", isofor, " ", input$ionmode, " ", input$adduct, wait=FALSE))
    
  })
  
  observeEvent(input$extractmzml, {
    mzFile <- as.character(input$entermzml)
    rt <- as.character(input$rtime)
    system(paste("python3 getSpectra.py ", mzFile, " ", rt, wait=FALSE))
  })
  
  observeEvent(input$calcperc, {
    system(paste("python3 CalcAbundancePerc.py compound_spectra.csv mz_calc.csv", wait=FALSE))
    
  })
  
})