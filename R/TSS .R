
#Plot TSSEnrichment

plotTSSEnrichment <- function(
  QCscATACRProj = NULL,
  groupBy = "Sample",
  chromSizes = getChromSizes(QCscATACRProj),
  TSS = getTSS(QCscATACRProj),
  flank = 2000,
  norm = 100,
  smooth = 11,
  pal = NULL,
  returnDF = FALSE,
  threads = getQCThreads(),
  logFile = createLogFile("plotTSSEnrichment")
){

  .validInput(input = QCscATACRProj, name = "QCscATACRProj", valid = c("QCscATACRProj"))
  .validInput(input = TSS, name = "TSS", valid = c("granges"))
  .validInput(input = flank, name = "flank", valid = c("integer"))
  .validInput(input = norm, name = "norm", valid = c("integer"))
  .validInput(input = smooth, name = "smooth", valid = c("integer"))
  .validInput(input = returnDF, name = "returnDF", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotTSSEnrichment Input-Parameters", logFile = logFile)

  chr <- paste0(seqnames(chromSizes))
  chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
  TSS <- sort(sortSeqlevels(TSS))
  splitTSS <- split(resize(TSS,1,"start"), seqnames(TSS))[chr]
  window <- 2 * flank + 1
  groups <- getCellColData(ArchRProj = QCscATACRProj, select = groupBy, drop = FALSE)
  uniqGroups <- gtools::mixedsort(unique(groups[,1]))

  if(threads > 1){
    h5disableFileLocking()
  }

  dfTSS <- .safelapply(seq_along(uniqGroups), function(x){

    .logDiffTime(paste0(uniqGroups[x], " Computing TSS (",x," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    cellx <- rownames(groups)[which(paste0(groups[,1]) == uniqGroups[x])]

    for(i in seq_along(chr)){

      TSSi <- splitTSS[[chr[i]]]

      covi <- unlist(suppressMessages(getFragmentsFromProject(
        QCscATACRProj = QCscATACRProj,
        subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% chr[i]],
        cellNames = cellx,
        logFile = logFile
      )), use.names=FALSE) %>%
        sort %>%
        {coverage(IRanges(c(start(.), end(.)), width = 1))}

      .logThis(covi, paste0(uniqGroups[x], " : Cov : ", chr[i]), logFile = logFile)

      if(i == 1){
        sumTSS <- rleSumsStranded(list(chr1=covi), list(chr1=TSSi), window, as.integer)
      }else{
        sumTSS <- sumTSS + rleSumsStranded(list(chr1=covi), list(chr1=TSSi), window, as.integer)
      }

      .logThis(sumTSS, paste0(uniqGroups[x], " : SumTSS : ", chr[i]), logFile = logFile)

    }

    normBy <- mean(sumTSS[c(1:norm,(flank*2-norm+1):(flank*2+1))])

    df <- DataFrame(
      group = uniqGroups[x],
      x = seq_along(sumTSS) - flank - 1,
      value = sumTSS,
      normValue = sumTSS / normBy,
      smoothValue = .centerRollMean(sumTSS/normBy, 11)
    )

    .logThis(df, paste0(uniqGroups[x], " : TSSDf"), logFile = logFile)

    .logDiffTime(paste0(uniqGroups[x], " Finished Computing TSS (",x," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    df

  }, threads = threads) %>% Reduce("rbind", .)

  .logThis(dfTSS, paste0("All : TSSDf"), logFile = logFile)

  .endLogging(logFile = logFile)

  if(threads > 1){
    h5enableFileLocking()
  }

  if(returnDF){

    return(dfTSS)

  }else{

    plotDF <- data.frame(x=dfTSS$x,v=dfTSS$smoothValue,group=dfTSS$group)
    plotDF <- plotDF[sort(unique(c(1,seq(1,nrow(plotDF),11),nrow(plotDF)))), , drop = FALSE]

    if(is.null(pal)){
      pal <- paletteDiscrete(values=unique(plotDF$group))
    }

    p <- ggplot(plotDF, aes(x,v,color=group)) +
      geom_line(size = 1) +
      theme_ArchR() +
      xlab("Distance From Center (bp)") +
      ylab("Normalized Insertion Profile") +
      scale_color_manual(values=pal) +
      scale_y_continuous(limits = c(0, max(plotDF$v)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0,0))

    p

  }

}

