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
library(dplyr)
shinyServer(function(input, output, session) {
    observeEvent(input$runimputation, {
      impute.df <- read.csv(input$psmfilename, header=TRUE)
      startAbundance <- as.integer(grep(paste(input$abundancecolumn,"$",sep=''), colnames(impute.df)))-1
      print(input$imputationlevel)
      system(paste("python3 imputationMV.py ", input$psmfilename, " ", input$replicatenum, " ", startAbundance, " ", input$imputationlevel, wait=FALSE))
      
    })
    dataFrame <- reactive({

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
      annotmir5a6.df <- xmir5a6.df[, c(as.integer(grep("Annotated.Sequence", colnames(xmir5a6.df))), as.integer(grep("Master.Protein.Accessions", colnames(xmir5a6.df))))]
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
      endAbundance <- as.integer(grep(paste(input$begchannels,"$",sep=''), colnames(xmir5a6.df)))+as.integer(input$numchannels)-1
      startAbundance <- as.integer(grep(paste(input$begchannels,"$",sep=''), colnames(xmir5a6.df)))
      isolint <- as.integer(grep("^Isolation.Interference", colnames(xmir5a6.df)))
      outfile <- input$outputfile
      pcaCtl <- input$pcacontrol
      pcaTreat <- input$pcatreatment
    
      #CHANGE df, peptix, fileIDix, isolinterfix, lessperc, startix, endix
      #isolinterfix, Isolation Interference [%] column
      #lessperc, is float for setting coisolation interference threshold (i.e. default 70.0)
      mir5a6.df <- prepareData_TS_PDPSM(xmir5a6.df, as.integer(grep("Annotated.Sequence", colnames(xmir5a6.df))), as.integer(grep("File.ID", colnames(xmir5a6.df))), isolint, as.numeric(input$lessperc), startAbundance, endAbundance)

      ymir5a6.lst <- separate_PDPSM(mir5a6.df, 2)

      xmir5a6.lst <- rmAnyMissing(ymir5a6.lst)
  
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
        
        return(mss)
      }
      resmir5a6.lst <- prepBlkAnnot(resmir5a6.df, 'mir5a6')
      resmir5a6.mss <- makeBlkMSS(resmir5a6.lst, 'mir5a6')
      # NORMALIZATION check with boxplot
     
      resmir5a6vsn.mss <- normalise(resmir5a6.mss, 'vsn')

      tiff(paste("Figures/NormalizationBoxPlot.tiff"), width = 4, height = 4, units = 'in', res=600)
      .plot(resmir5a6vsn.mss)
      dev.off()
      
      pd <- phenoData(resmir5a6vsn.mss)$TreatmentGroup
      names(pd) <- sampleNames(resmir5a6vsn.mss)
      
      
      plotPCA_sc_v2 <- function(m, pdat, component, title='') { # select components
        ## component: either 1 (comp1vscomp2) or 2 (comp2vscomp3)
        pca <- prcomp(m, scale=FALSE)
        df <- as.data.frame(pca$rotation[, 1:4])
        df <- namerows(df, col.name='Samples')
        
        spl <- df$Samples
        cl <- pdat[match(spl, names(pdat))]
        spl <- ifelse(cl==1, pcaCtl, pcaTreat)
        df$Samples <- spl
        
        if (component=='1') { 
          p <- ggplot(df, aes(PC1, PC2)) + geom_point(shape=21, size = 2, stroke=1, color="black", aes(fill=Samples)) + scale_fill_manual(values=c("white", "black"))
        } else if (component=='2') {
          p <- ggplot(df, aes(PC2, PC3, colour=Samples)) + geom_point(size=2) 
          
        }
        
        p <- p + theme(legend.position='right', legend.title=element_blank())
        p <- p + labs(title=title)
        
        return(p)
      }
      
      #change file name
      e <- exprs(resmir5a6vsn.mss)
      p <- plotPCA_sc_v2(e, pd, '1', title=paste('', '')) +
        theme_classic()
      tiff(paste("Figures/PCAplot.tiff"), width = 4, height = 4, units = 'in', res = 600)
      
      plot(p)
      dev.off()
      
      
      group <- factor(phenoData(resmir5a6vsn.mss)$TreatmentGroup)
      
      design <- model.matrix(~0+group)
      colnames(design) <- c('Ctrl', 'Transgn')

      fit <- lmFit(e, design)
      cm <- makeContrasts(Ctrl-Transgn, levels=design)
      
      fit2 <- contrasts.fit(fit, cm)
      fit2 <- eBayes(fit2)

      ttUp.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
      ttUp.df$symbol <- unlist(mget(rownames(ttUp.df), uniprotmir5a62sym, ifnotfound=rownames(ttUp.df)))
      ttUp.df$FC <- ifelse(ttUp.df$logFC >= 0, inv.glog2(ttUp.df$logFC), -inv.glog2(-ttUp.df$logFC))
      ttUp.df$logPval <- -log10(ttUp.df[,c(2)])
      ttUp.df <- ttUp.df[which(ttUp.df$logFC  >= 0.58 & ttUp.df$logPval >= 1.3),]
      write.table(ttUp.df, file=paste("Figures/StatsUpregulated.csv"), quote=FALSE, sep=',', row.names = FALSE)
      rm(ttUp.df)

      ttDown.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
      ttDown.df$symbol <- unlist(mget(rownames(ttDown.df), uniprotmir5a62sym, ifnotfound=rownames(ttDown.df)))
      ttDown.df$FC <- ifelse(ttDown.df$logFC >= 0, inv.glog2(ttDown.df$logFC), -inv.glog2(-ttDown.df$logFC))
      ttDown.df$logPval <- -log10(ttDown.df[,c(2)])
      ttDown.df <- ttDown.df[which(ttDown.df$logFC  <= -0.58 & ttDown.df$logPval >= 1.3),]
      write.table(ttDown.df, file=paste("Figures/StatsDownregulated.csv"), quote=FALSE, sep=',', row.names = FALSE)
      rm(ttDown.df)

      tt.df <- topTable(fit2, number=Inf, sort.by ='p', p.value=1)[, c(1, 4, 5)]
      tt.df$symbol <- unlist(mget(rownames(tt.df), uniprotmir5a62sym, ifnotfound=rownames(tt.df)))
      tt.df$FC <- ifelse(tt.df$logFC >= 0, inv.glog2(tt.df$logFC), -inv.glog2(-tt.df$logFC))
      tt.df$logPval <- -log10(tt.df[,c(2)])
      
      write.table(tt.df, file=paste("Figures/StatsTable.csv"), quote=FALSE, sep=',', row.names = FALSE)
      
        tt.df
        
      
    })
    
    dataFilter <- reactive({
      dataFrame()[dataFrame()$logFC > input$fcCut[1] & dataFrame()$logFC < input$fcCut[2],]
    })
    
    output$volcanoPlot <- renderPlot({
        highlight_df <- dataFilter() %>% 
          filter(symbol==input$protint)
        highlight_df_down <- dataFilter() %>% 
          filter(logFC<=-0.58)
        highlight_df_up <- dataFilter() %>% 
          filter(logFC>=0.58)
        ggplot(dataFilter(),aes(x=logFC,y=logPval)) + geom_point(size=2, alpha=1, col='black') +
        labs(title=input$plottitle, x=input$xaxis, y=input$yaxis) +
        theme_update(plot.title=element_text(hjust=0.5), legend.position='none') +
        geom_point(data=dataFilter(), stat='identity', aes(colour=cut(logFC, c(-Inf,-0.58,0.58,5))), size=1) + geom_hline(yintercept=-log(0.05,10), linetype="3313", colour="grey") + geom_vline(xintercept=0.58, linetype="3313", colour="grey") + geom_vline(xintercept=-0.58, linetype="3313", colour="grey") +
          scale_color_manual(name = "logFC",
                             values = c("(-Inf,-0.58]" = "blue",
                                        "(-0.58,0.58]" = "gray",
                                        "(0.58,5]" = "red"),
                             labels = c("decreased", "insignificant", "increased")) +
        geom_point(data=highlight_df, aes(x=logFC,y=logPval), color='green',size=2,alpha=1, col='black') +
        #geom_text_repel(data=highlight_df, aes(x=logFC, y=logPval, label=highlight_df$symbol), colour='forestgreen', size=2) +
        geom_text_repel(data=highlight_df_down, aes(x=logFC, y=logPval, label=highlight_df_down$symbol), colour='black', size=2) +
        geom_text_repel(data=highlight_df_up, aes(x=logFC, y=logPval, label=highlight_df_up$symbol), colour='black', size=2) +
        theme_classic()
        

    })
    
    plotOutput <- reactive({
      print(input$protint)
      highlight_df <- dataFilter() %>% 
        filter (symbol==input$protint)
      highlight_df_down <- dataFilter() %>% 
        filter(logFC<=-0.58)
      highlight_df_up <- dataFilter() %>% 
        filter(logFC>=0.58)
      ggplot(dataFilter(),aes(x=logFC,y=logPval)) + geom_point(size=2, alpha=1, col='black') +
        labs(title=input$plottitle, x=input$xaxis, y=input$yaxis) +
        theme_update(plot.title=element_text(hjust=0.5), legend.position='none') +
        geom_point(data=dataFilter(), stat='identity', aes(colour=cut(logFC, c(-Inf,-0.58,0.58,5))), size=1) + geom_hline(yintercept=-log(0.05,10), linetype="3313", colour="grey") + geom_vline(xintercept=0.58, linetype="3313", colour="grey") + geom_vline(xintercept=-0.58, linetype="3313", colour="grey") +
        scale_color_manual(name = "logFC",
                           values = c("(-Inf,-0.58]" = "blue",
                                      "(-0.58,0.58]" = "gray",
                                      "(0.58,5]" = "red"),
                           labels = c("decreased", "insignificant", "increased")) +
        geom_point(data=highlight_df, aes(x=logFC,y=logPval,label=symbol), color='green',size=2, alpha=1, col='black') +
        #geom_text_repel(data=highlight_df, aes(x=logFC, y=logPval, label=highlight_df$symbol), colour='forestgreen', size=2) +
        geom_text_repel(data=highlight_df_down, aes(x=logFC, y=logPval, label=highlight_df_down$symbol), colour='black', size=2) +
        geom_text_repel(data=highlight_df_up, aes(x=logFC, y=logPval, label=highlight_df_up$symbol), colour='black', size=2) +
        theme_classic()
      
      
    })


    observeEvent(input$downloadPlot, {
        ggsave(paste("Figures/VolcanoPlot.png"),plotOutput(), width = 8, height = 4, dpi=600)
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

