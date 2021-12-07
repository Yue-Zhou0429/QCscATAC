#' Plot the fragment size distribution

plotFragmentSizes <- function(
  QCscATACRProj = NULL,
  groupBy = "Sample",
  chromSizes = getChromSizes(QCscATACRProj),
  maxSize = 750,
  pal = NULL,
  returnDF = FALSE,
  threads = getQCThreads(),
  logFile = createLogFile("plotFragmentSizes")
){

  .validInput(input = QCscATACRProj, name = "QCscATACRProj", valid = c("QCscATACRProj"))
  .validInput(input = maxSize, name = "maxSize", valid = c("integer"))
  .validInput(input = returnDF, name = "returnDF", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotFragmentSizes Input-Parameters", logFile = logFile)

  chr <- paste0(seqnames(chromSizes))
  groups <- getCellColData(QCscATACRProj = QCscATACRProj, select = groupBy, drop = FALSE)
  uniqGroups <- gtools::mixedsort(unique(groups[,1]))

  if(threads > 1){
    h5disableFileLocking()
  }

  dfFS <- .safelapply(seq_along(uniqGroups), function(x){

    .logDiffTime(paste0(uniqGroups[x], " Computing FragmentSizes (",x," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    cellx <- rownames(groups)[which(paste0(groups[,1]) == uniqGroups[x])]

    for(i in seq_along(chr)){
      if(i == 1){
        fsi <- unlist(suppressMessages(getFragmentsFromProject(
          ArchRProj = ArchRProj,
          subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% chr[i]],
          cellNames = cellx,
          logFile = logFile
        )), use.names=FALSE) %>% width %>% tabulate(nbins = maxSize)
      }else{
        fsi <- fsi + unlist(suppressMessages(getFragmentsFromProject(
          ArchRProj = ArchRProj,
          subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% chr[i]],
          cellNames = cellx,
          logFile = logFile
        )), use.names=FALSE) %>% width %>% tabulate(nbins = maxSize)
      }
      .logThis(fsi, paste0(uniqGroups[x], " : FragSizes : ", chr[i]), logFile = logFile)
    }

    df <- DataFrame(
      group = uniqGroups[x],
      fragmentSize = seq_along(fsi),
      fragmentPercent = round(100*fsi/sum(fsi),4)
    )

    .logThis(df, paste0(uniqGroups[x], " : Frag DF"), logFile = logFile)

    .logDiffTime(paste0(uniqGroups[x], " Finished Computing FragmentSizes (",x," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    df

  }, threads = threads) %>% Reduce("rbind", .)

  .logThis(dfFS, paste0("All : FragSizes DF"), logFile = logFile)

  .endLogging(logFile = logFile)

  if(threads > 1){
    h5enableFileLocking()
  }

  if(returnDF){

    return(dfFS)

  }else{

    plotDF <- data.frame(dfFS)

    if(is.null(pal)){
      pal <- paletteDiscrete(values=unique(plotDF$group))
    }

    p <- ggplot(plotDF, aes(fragmentSize, fragmentPercent,color=group)) +
      geom_line(size = 1) +
      theme_ArchR() +
      xlab("ATAC-seq Fragment Size (bp)") +
      ylab("Percentage of Fragments") +
      scale_color_manual(values=pal) +
      scale_y_continuous(limits = c(0, max(plotDF$fragmentPercent)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$fragmentSize), max(plotDF$fragmentSize)), expand = c(0,0))

    p

  }

}


