library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("superhero"),
    fileInput('csvfile', 'Input File'),
    fileInput('uniprotout', 'Uniprot File'),
    fileInput('unitogene', 'Uniprot and Gene Name File'),
    textInput('peptideCol', 'Peptide Column', value = 5),
    textInput('accessionCol', 'Accession Column', value = 8),
    textInput('fileIDix', 'File ID Column', value = 29),
    textInput('isolinterfix', 'Isolation Interference Column', value = 23),
    textInput('lessperc', 'Coisolation Interference Threshold (default 70%)', value = 70.0),
    textInput('startix', 'Abundance Start Column', value = 30),
    textInput('endix', 'Abundance End Column', value = 39),
    radioButtons('channel126', 'Channel 126', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel127N', 'Channel 127N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel127C', 'Channel 127C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel128N', 'Channel 128N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel128C', 'Channel 128C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel129N', 'Channel 129N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel129C', 'Channel 129C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel130N', 'Channel 130N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel130C', 'Channel 130C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel131N', 'Channel 131N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel131C', 'Channel 131C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel132N', 'Channel 132N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel132C', 'Channel 132C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel133N', 'Channel 133N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel133C', 'Channel 133C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel134N', 'Channel 134N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    radioButtons('channel134C', 'Channel 134C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "Treatment 2",
      "NA"
    ),
    choiceValues = list(
      '0','1','2','3'
    )),
    textInput('outputfile', 'Output File Name', value = 'Control_vs_Treatment'),
    actionButton('buttonId', 'run script'),
    titlePanel("Volcano Plot"),
    plotOutput('volcanoPlot',click='plot_click'),
    sliderInput('fcCut', label="log(FC) cutoff",min=-2,max=2,value=c(-2,-2), step=0.1, width="600px"),

    #here the table for the clicked points:
    tableOutput('clickedPoints')
  )
  
)