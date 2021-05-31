library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("cyborg"),
    
                  textInput('entermzml', 'Enter mzML File'),
                  textInput('enterpsm', 'Enter PSM File'),
                  actionButton('extractmzml', 'Extract mzML.'),
                  actionButton('calcperc', 'Extract Peak Areas.'),
))