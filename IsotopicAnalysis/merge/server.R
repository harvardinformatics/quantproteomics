options(shiny.maxRequestSize=300*1024^2) 

library(mzR)
library(mwshiny)

shinyServer(function(input, output, session) {
  observeEvent(input$calciso, {
    isofor <- as.character(input$isoformula)
    print(isofor)
    system(paste("python3 IsotopeAnalysis_Merge.py ", isofor, " ", input$ionmode, " ", input$adduct, " ", mzFile, " ", wait=FALSE))
    
  })
  

  
})