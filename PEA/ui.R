library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("superhero"),
    textInput('psmfilename', 'PSM Filename'),
    textInput('replicatenum1', 'Replicate Number of Sample 1'),
    textInput('replicatenum2', 'Replicate Number of Sample 2'),
    textInput('abundancecolumn', 'Abundance Column'),
    radioButtons('imputationlevel', 'Level of Imputation', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "High",
      "Medium",
      "Low"
    ),
    choiceValues = list(
      '2','1','0'
    )),
    actionButton('runimputation', 'Impute my Missing Values!'),
    textInput('PSMfile', 'PSM Imputed File'),
    textInput('Protfile', 'PD Protein File'),
    actionButton('runPDfilter', 'Use PD filters'),
    fileInput('csvfile', 'Input File'),
    fileInput('uniprotout', 'Uniprot File'),
    fileInput('unitogene', 'Uniprot and Gene Name File'),
    textInput('protnorm', 'Protein to Normalize to', value = 'NA'),
    textInput('lessperc', 'Coisolation Interference Threshold (default 70%)', value = 70.0),
    textInput('plottitle', 'Plot Title', value = 'Differential Expressed Proteins for Treatment/Control at P Value <= 0.05'),
    textInput('xaxis', 'Plot X-axis', value = 'log2(Treatment/Control)'),
    textInput('yaxis', 'Plot Y-axis', value = '-log10(nominalpval)'),
    textInput('pcacontrol', 'PCA Control', value = 'Control'),
    textInput('pcatreatment', 'PCA Treatment', value = 'Treatment'),
    textInput('protint1', 'Protein of Interest 1', value = 'NA'),
    textInput('protint2', 'Protein of Interest 2', value = 'NA'),
    textInput('protint3', 'Protein of Interest 3', value = 'NA'),
    radioButtons('channel126', 'Channel 126', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel127N', 'Channel 127N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel127C', 'Channel 127C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel128N', 'Channel 128N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel128C', 'Channel 128C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel129N', 'Channel 129N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel129C', 'Channel 129C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel130N', 'Channel 130N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel130C', 'Channel 130C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel131N', 'Channel 131N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel131C', 'Channel 131C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel132N', 'Channel 132N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel132C', 'Channel 132C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel133N', 'Channel 133N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel133C', 'Channel 133C', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),
    radioButtons('channel134N', 'Channel 134N', inline = TRUE, choices = NULL, selected = NULL, choiceNames = list(
      "Control",
      "Treatment",
      "NA"
    ),
    choiceValues = list(
      '1','0','2'
    )),

    
    actionButton('buttonId', 'run script'),
    titlePanel("Volcano Plot"),
    plotOutput('volcanoPlot',click='plot_click'),
    sliderInput('fcCut', label="log(FC) cutoff",min=-3,max=3,value=c(-3,3), step=0.1, width="600px"),
    actionButton('downloadPlot', 'Download Plot'),
    #here the table for the clicked points:
    tableOutput('clickedPoints')
  )
  
)