normalise_MSnSet <- function(object, method, ...) {
  if (method == "vsn") {
    e <- exprs(vsn2(exprs(object), ...))
  } else if (method == "quantiles") {
    e <- preprocessCore::normalize.quantiles(exprs(object), ...)
  } else if (method == "quantiles.robust") {
    e <- preprocessCore::normalize.quantiles.robust(exprs(object), ...)
  } else if (method == "center.mean") {
    e <- exprs(object)
    center <- colMeans(e, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE, ...)
  } else if (method == "center.median") {
    e <- exprs(object)
    center <- apply(e, 2L, median, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE, ...)
  } else if (method == "diff.median") {
      e <- exprs(object)
      med <- median(as.numeric(e), na.rm = TRUE)
      cmeds <- apply(e, 2L, median, na.rm = TRUE)
      e <- sweep(e, 2L, cmeds - med)
  } else {
    switch(method,
           max = div <- .rowMaxs(exprs(object), na.rm = TRUE),
           sum = div <- rowSums(exprs(object), na.rm = TRUE))
    e <- exprs(object)/div
  }
  rownames(e) <- rownames(exprs(object))
  colnames(e) <- colnames(exprs(object))
    exprs(object) <- e
  object@processingData@processing <- c(object@processingData@processing, paste("Normalised (", method ,"): ", date(), sep = ""))
  object@processingData@normalised <- TRUE
  if (validObject(object))
    return(object)
}

.plot <- function(x,ttl=NULL) {
    boxplot(exprs(x),
          main=ifelse(is.null(ttl),processingData(x)@processing[2],ttl),
          cex.main=.8,
          cex.lab=0.7,
          cex.axis=.7, # 0.9 for axis label
          cex=0.7, las=2, pch=20)
    grid()
}

plotPCA_sc_v2 <- function(m, pdat, component, title='') { # select components
    ## component: either 1 (comp1vscomp2) or 2 (comp2vscomp3)
    pca <- prcomp(m, scale=FALSE)
    df <- as.data.frame(pca$rotation[, 1:4])
    df <- namerows(df, col.name='Samples')
    
    spl <- df$Samples
    cl <- pdat[match(spl, names(pdat))]
    spl <- ifelse(cl==0, 'Ctrl', 'Transg')
    df$Samples <- spl

    if (component=='1') { 
        p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2) + scale_color_grey()
    } else if (component=='2') {
        p <- ggplot(df, aes(PC2, PC3, colour=Samples)) + geom_point(size=2) + scale_color_grey()
        #p <- ggplot(df, aes(PC3, PC4, colour=Samples)) + geom_point(size=2)
    }
    
    p <- p + theme(legend.position='right', legend.title=element_blank())
    p <- p + labs(title=title)
    
    return(p)
}

prepareData_TS_PDPSM <- function(df, peptix, fileIDix, isolinterfix, lessperc, abu) { # TS: Targeted and SCOPED
    # lessperc: float for setting coisolation interference threshold
    filterIsolationInterference <- function(adf, ix, perc) {
        # D040219 - remove rows with missing interference percentage
        perc <- as.numeric(perc)
        ix <- as.integer(ix)
        isolinterf <- adf[, ix]
        adf <- adf[isolinterf < perc & !is.na(isolinterf), ]
        return(adf)
    }
    filterCarrierEmptyChannels_v0 <- function(adf, carrier, empty) {
        cix <- grep(paste(carrier, empty, sep='|'), colnames(adf))
        adf <- adf[, -cix]
        return(adf)
    }
    filterCarrierEmptyChannels <- function(adf, carrierORempty, ...) {
        alst <- list(...)
        carrierORempty <- c(carrierORempty, unlist(alst))
        cix <- grep(paste(carrierORempty, collapse='|'), colnames(adf))
        adf <- adf[, -cix]
        return(adf)
    }
    rmEmptyRows <- function(df) {
        emptyrows <- apply(df, 1, function(x) all(is.na(x)))
        
        return(df[!emptyrows, ])
    }
    cleanColumnNames <- function(df) {
        cn  <-  colnames(df)
        cn <- sub('Abundance\\.\\.', '_', colnames(df))
        colnames(df) <- cn
        
        return(df)
    }
    rmPeptInMultipleProt  <- function(df) {
        return(df[df$'X..Proteins' == 1, ])
    }
    df <- cleanColumnNames(df) 
    df <- filterIsolationInterference(df, isolinterfix, lessperc)
    df <- rmPeptInMultipleProt(df)

    peptix <- as.integer(peptix)
    fileIDix <- as.integer(fileIDix)
    #startix <- as.integer(startix)
    #endix <- as.integer(endix)
    df <- df[, c(peptix, fileIDix, abu)]
    
    # CAREFUL - this should be moved into the argument list!
    #df <- filterCarrierEmptyChannels(df, '126', '127')
    
    colnames(df)[1] <- 'PepSeq'
    
    return(df)
}

separate_PDPSM <- function(df, runix) {
    runix <- as.integer(runix)
    splt.lst <- split(df, factor(df[, runix, drop=TRUE]))
    
    splt.lst <- lapply(splt.lst, function(x) {
        prefix <- unique(x[, runix, drop=TRUE])
        abundix <- grep('[[:digit:]]', colnames(x))
        colnames(x)[abundix] <- paste(prefix, colnames(x)[abundix], sep='')
        colnames(x)[abundix] <- sub('Abundance\\.\\.', '', colnames(x)[abundix])
        x <- x[-runix]
    })
    
    return(splt.lst)
}

rmAnyMissing <- function(lst) {
    checked.lst <- lapply(lst, function(df) {
        discard <- apply(df, 1, function(x) any(is.na(x[-1])))
        return(df[!discard, ])
    })
    
    return(checked.lst)
}
