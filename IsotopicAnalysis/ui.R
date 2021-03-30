library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("cyborg"),
                  textInput('isoformula', 'Isotope Formula'),
                  textInput('ionmode', 'Ion Mode'),
                  textInput('adduct', 'Adduct'),
                  actionButton('calciso', 'Calculate Isotope.'),
                  textInput('entermzml', 'Enter mzML File'),
                  actionButton('extractmzml', 'Extract mzML.'),
))