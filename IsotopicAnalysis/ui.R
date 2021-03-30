library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("superhero"),
                  textInput('isoformula', 'Isotope Formula'),
                  textInput('ionmode', 'Ion Mode'),
                  textInput('adduct', 'Adduct'),
                  actionButton('calciso', 'Calculate Isotope.'),
))