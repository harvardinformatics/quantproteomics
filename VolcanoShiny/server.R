options(shiny.maxRequestSize=300*1024^2) 

library(MSnbase)
library(reshape)
library(lattice)
library(ggplot2)
library(limma)

library(ROTS)

library(seqinr)

library(ggrepel)
library(MKmisc) # glog2

library(gplots)


library(gtools)
library(magick)

library(mwshiny)
shinyServer(function(input, output, session) {

    dataFrame <- reactive({
      X11.options(width=5, height=5, xpos=1200, ypos=500)
      options(editor="/usr/bin/vim")
      options(stringsAsFactors=FALSE)
      
      source('functions_for_proteomics_Rcode.R')
      
      # LOAD DATA
      #CHANGE input file name
      inFile<-input$csvfile
      if (is.null(inFile))
        return(NULL)
      xmir5a6.df <- read.csv(inFile$datapath, header=TRUE)
    
      print('printing args from R code') 
      if (!exists("args")) {
        suppressPackageStartupMessages(library("argparse"))
        parser <- ArgumentParser()
        parser$add_argument("-a", "--arg1", type="character", defalt="a",
                            help="First parameter [default %(defult)s]")
        parser$add_argument("-b", "--arg2", type="character", defalt="b",
                            help="Second parameter [default %(defult)s]")
        args <- parser$parse_args()
      }
      
      
      # annotation
      annotmir5a6.df <- xmir5a6.df[, c(as.integer(input$peptideCol), as.integer(input$accessionCol))]
      colnames(annotmir5a6.df) <- c('Pept', 'Acc')
      annotmir5a6.df$Acc <- as.character(sapply(annotmir5a6.df$Acc, function(x) unlist(strsplit(x, split=';'))[1]))
      annotmir5a6.df <- annotmir5a6.df[!is.na(annotmir5a6.df$Acc), ]
      
      # lookup
      pepseqmir5a62acc <- new.env(hash=TRUE)
      apply(annotmir5a6.df, 1, function(x) {
        pepseqmir5a62acc[[x[1]]] <- x[2]
      })
      
      uprot1<-input$uniprotout
      if (is.null(uprot1))
        return(NULL)
      uprot2<-input$unitogene
      if (is.null(uprot2))
        return(NULL)
      write.table(unique(annotmir5a6.df$Acc), uprot1$datapath, quote=FALSE, sep=',', row.names=FALSE, col.names=FALSE)
      # download annotation from UniProt
      uniprot2genename.df <- read.table(uprot2$datapath, header=FALSE, sep=',', quote='')
      uniprotmir5a62sym <- new.env(hash=TRUE)
      apply(uniprot2genename.df, 1, function(x) {
        x <- as.character(x)
        uniprotmir5a62sym[[x[1]]] <- x[2]
      })
      
      #CHANGE df, peptix, fileIDix, isolinterfix, lessperc, startix, endix
      #isolinterfix, Isolation Interference [%] column
      #lessperc, is float for setting coisolation interference threshold (i.e. default 70.0)
      mir5a6.df <- prepareData_TS_PDPSM(xmir5a6.df, as.integer(input$peptideCol), as.integer(input$fileIDix), as.integer(input$isolinterfix), as.numeric(input$lessperc), as.integer(input$startix), as.integer(input$endix))
      ymir5a6.lst <- separate_PDPSM(mir5a6.df, 2)
      xmir5a6.lst <- rmAnyMissing(ymir5a6.lst)
      
      #xmir5a6.lst <- xmir5a6.lst[paste('F1', seq(20), sep='.')]
      
      # add protein column
      mir5a6.lst <- lapply(xmir5a6.lst, function(df) {
        rownames(df) <- make.unique(df$PepSeq, sep=';')    
        df$Prot <- as.character(unlist(mget(df$PepSeq, pepseqmir5a62acc, ifnotfound=NA)))
        df <- df[!is.na(df$Prot), ]
        df <- df[-1]
        df <- df[, c(ncol(df), seq(ncol(df)-1))]
        o <- order(df$Prot)
        df <- df[o, ]
        return(df)
      })
      
      mir5a6.lst <- lapply(mir5a6.lst, function(df) {
        colnames(df)  <- sub('F.*_', 'F_', colnames(df))
        return(df)
      })
      
      xresmir5a6.df <- do.call(rbind, mir5a6.lst)
      resmir5a6.df <- aggregate(xresmir5a6.df[-1], xresmir5a6.df[1], sum) # aggregating! maybe not?
      
      # MAKE MSnbase OBJECT
      prepBlkAnnot <- function(df, suff) {
        adf <- df
        adf <- adf[order(adf$Prot, decreasing=FALSE), ]
        
        fdf <- data.frame(ID=adf$Prot, Acc=adf$Prot)
        rownames(fdf) <- fdf$ID
        
        rownames(adf) <- adf$Prot
        bm <- adf[-1]
        bm[colnames(bm)] <- sapply(bm[colnames(bm)], as.numeric)
        bm <- as.matrix(bm)
        
        bm.cnames <- sapply(colnames(bm), function(x) {
          if (grepl('126', x)) {
            x <- paste(x, as.integer(input$channel126), sep=',')
          } else if (grepl('127N', x)) {
            x <- paste(x, as.integer(input$channel127N), sep=',')
          } else if (grepl('127C', x)) {
            x <- paste(x, as.integer(input$channel127C), sep=',')
          } else if (grepl('128N', x)) {
            x <- paste(x, as.integer(input$channel128N), sep=',')
          } else if (grepl('128C', x)) {
            x <- paste(x, as.integer(input$channel128C), sep=',')
          } else if (grepl('129N', x)) {
            x <- paste(x, as.integer(input$channel129N), sep=',')
          } else if (grepl('129C', x)) {
            x <- paste(x, as.integer(input$channel129C), sep=',')
          } else if (grepl('130N', x)) {
            x <- paste(x, as.integer(input$channel130N), sep=',')
          } else if (grepl('130C', x)) {
            x <- paste(x, as.integer(input$channel130C), sep=',')
          } else if (grepl('131N', x)) {
            x <- paste(x, as.integer(input$channel131N), sep=',')
          } else if (grepl('131C', x)) {
            x <- paste(x, as.integer(input$channel131C), sep=',')
          } else if (grepl('132N', x)) {
            x <- paste(x, as.integer(input$channel132N), sep=',')
          } else if (grepl('132C', x)) {
            x <- paste(x, as.integer(input$channel132C), sep=',')
          } else if (grepl('133N', x)) {
            x <- paste(x, as.integer(input$channel133N), sep=',')
          } else if (grepl('133C', x)) {
            x <- paste(x, as.integer(input$channel133C), sep=',')
          } else if (grepl('134N', x)) {
            x <- paste(x, as.integer(input$channel134N), sep=',')
          } else if (grepl('134C', x)) {
            x <- paste(x, as.integer(input$channel134C), sep=',')
          } 
        })
        write.table(bm.cnames, file=paste('pData_', suff, '.txt', sep=''), col.names=paste('TreatmentGroup', sep=','),
                    row.names=FALSE, quote=FALSE)
        
        return(list(bm, fdf))
      }
      
      makeBlkMSS <- function(lst, suff) {
        pd <- read.csv(paste('pData_', suff, '.txt', sep=''))
        mss <- MSnSet(lst[[1]], lst[[2]], pd)
        #mss <- mss[, grep('126|127', sampleNames(mss), invert=TRUE)]
        
        return(mss)
      }
      resmir5a6.lst <- prepBlkAnnot(resmir5a6.df, 'mir5a6')
      resmir5a6.mss <- makeBlkMSS(resmir5a6.lst, 'mir5a6')
      
      # NORMALIZATION check with boxplot
      #change file name
      resmir5a6vsn.mss <- normalise(resmir5a6.mss, 'vsn')
      tiff(paste("/Users/ald533/Desktop/ProductionAndInformatics/FASInformatics/Figures/", input$outputfile, "_NormalizationBox.tiff"), width = 4, height = 4, units = 'in', res=600)
      .plot(resmir5a6vsn.mss)
      dev.off()
      
      pd <- phenoData(resmir5a6vsn.mss)$TreatmentGroup
      names(pd) <- sampleNames(resmir5a6vsn.mss)
      
      #change file name
      e <- exprs(resmir5a6vsn.mss)
      p <- plotPCA_sc_v2(e, pd, '1', title=paste('', '')) +
        theme_classic()
      tiff(paste("/Users/ald533/Desktop/ProductionAndInformatics/FASInformatics/Figures/", input$outputfile, "_PCA.tiff"), width = 4, height = 4, units = 'in', res = 600)
      
      plot(p)
      dev.off()
      
      
      group <- factor(phenoData(resmir5a6vsn.mss)$TreatmentGroup)
      design <- model.matrix(~0+group)
      colnames(design) <- c('Ctrl', 'Transgn')
      
      fit <- lmFit(e, design)
      
      cm <- makeContrasts(Ctrl-Transgn, levels=design)
      fit2 <- contrasts.fit(fit, cm)
      fit2 <- eBayes(fit2)
      
  
      tt.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
      tt.df$symbol <- unlist(mget(rownames(tt.df), uniprotmir5a62sym, ifnotfound=rownames(tt.df)))
      tt.df$FC <- ifelse(tt.df$logFC >= 0, inv.glog2(tt.df$logFC), -inv.glog2(-tt.df$logFC))
      tt.df$logPval <- -log10(tt.df[,c(2)])
      write.table(tt.df, file=paste("/Users/ald533/Desktop/ProductionAndInformatics/FASInformatics/Figures/", input$outputfile, '_Stats.csv'), quote=FALSE, sep='\t')
      
      
        tt.df
    })
    
    dataFilter <- reactive({
      dataFrame()[dataFrame()$logFC > input$fcCut[1] & dataFrame()$logFC < input$fcCut[2],]
    })
    
    output$volcanoPlot <- renderPlot({
        
        ggplot(dataFilter(),aes(x=logFC,y=logPval)) + geom_point(size=2, alpha=1, col='black') +
        labs(title='Diff Expressed Proteins for Control vs Treatment at P Value <= 0.05', x='log(FC)', y='-log(nominalpval)') +
        theme(plot.title=element_text(size=10, vjust=1, hjust=0.5), legend.position='none') +
        geom_point(data=dataFilter(), stat='identity', aes(colour=cut(logPval, c(0,1.3,Inf))), size=1) + geom_hline(yintercept=-log(0.05,10)) +
          scale_color_manual(name = "-log(P-value)",
                             values = c("(0,1.3]" = "blue",
                                        "(1.3,Inf]" = "red"),
                             labels = c("insignificant", "significant")) +
        geom_text_repel(data=dataFilter()[1:30, ], aes(x=logFC, y=logPval, label=dataFilter()$symbol[1:30]), colour='forestgreen', size=2) +
        theme_classic()
        

    })
    
    plotOutput <- reactive({
      
      ggplot(dataFilter(),aes(x=logFC,y=logPval)) + geom_point(size=2, alpha=1, col='black') +
        labs(title='Diff Expressed Proteins for Control vs Treatment at P Value <= 0.05', x='log(FC)', y='-log(nominalpval)') +
        theme(plot.title=element_text(size=10, vjust=1, hjust=0.5), legend.position='none') +
        geom_point(data=dataFilter(), stat='identity', aes(colour=cut(logPval, c(0,1.3,Inf))), size=1) + geom_hline(yintercept=-log(0.05,10)) +
        scale_color_manual(name = "-log(P-value)",
                           values = c("(0,1.3]" = "blue",
                                      "(1.3,Inf]" = "red"),
                           labels = c("insignificant", "significant")) +
        geom_text_repel(data=dataFilter()[1:30, ], aes(x=logFC, y=logPval, label=dataFilter()$symbol[1:30]), colour='forestgreen', size=2) +
        theme_classic()
      
      
    })


    observeEvent(input$downloadPlot, {
        ggsave(paste("/Users/ald533/Desktop/ProductionAndInformatics/FASInformatics/Figures/", input$outputfile, '_Volcano.png'),plotOutput())
      })
    
    

    clicked <- reactive({
      # We need to tell it what the x and y variables are:
      nearPoints(dataFilter(), input$plot_click, xvar = "logFC", yvar = "logPval")
    })
    
    #output those points into a table
    output$clickedPoints <- renderTable({
      clicked()
    }, rownames = T)
    
    


})

