setwd("/panfs/pfs.local/work/colella/ben/002.simulations/003.SecondaryContact_vs_IBD/03.conStruct/")

library(ggplot2) 
library(conStruct)
library(scales)
library(reshape2)

#############################################################################
# Define custom structure plot functions, based on structurePlot from strataG
#############################################################################


betterStructurePlot <- function(
  q.mat, 
  ind.col = 1,
  pop.col = 3, 
  prob.col = 4, 
  sort.probs = TRUE,
  label.pops = FALSE,
  label.inds = TRUE,
  col = NULL, 
  horiz = TRUE, 
  type = NULL,
  legend.position = c("top", "left", "right", "bottom", "none"),
  plot = TRUE
) {
  
  legend.position <- match.arg(legend.position)
  
  # convert q.mat to sorted data.table
  prob.cols <- prob.col:ncol(q.mat)
  qm <- as.data.frame(q.mat)[, c(ind.col, pop.col, prob.cols), drop = FALSE]
  qm[, 2] <- factor(
    qm[, 2], 
    levels = sort(unique(qm[, 2]), decreasing = !horiz)
  )
  sort.cols <- c(2, if(sort.probs) 3:ncol(qm) else NULL)
  i <- do.call(
    order, 
    c(as.list(qm[, sort.cols, drop = FALSE]), decreasing = TRUE)
  )
  qm <- qm[i, ] # this reverses the order of pops for plotting
  qm$x <- 1:nrow(qm)
  
  # Get population frequencies, centers and dividing points
  pop.freq <- table(qm[, 2])
  levels(qm[, 2]) <- paste(
    levels(qm[, 2]), "\n(n = ", pop.freq, ")", sep = ""
  )
  pop.cntr <- tapply(qm$x, qm[, 2], mean)
  ind.cntr <- tapply(qm$x, qm[, 1], mean)
  pop.div <- rev(tapply(qm$x, qm[, 2], min))[-1] - 0.5
  
  # Create data.frame for plotting
  df <- melt(qm[, c("x", colnames(qm)[-ncol(qm)])], id.vars = c(1:3),
             variable.name = "Group", value.name = "probability")
  colnames(df)[1:3] <- c("x", "id", "population")
  df <- df[order(-as.numeric(df$Group), df$probability), ]
  
  type <- if(is.null(type)) {
    if(nrow(df) <= 100) "bar" else "area"
  } else {
    match.arg(type, c("bar", "area"))
  }
  
  # Plot stacked bar graphs
  g <- ggplot2::ggplot(df, ggplot2::aes_string("x", "probability")) +  
    switch(
      type,
      area = ggplot2::geom_area(
        ggplot2::aes_string(fill = "Group"), 
        stat = "identity"
      ),
      bar = ggplot2::geom_bar(
        ggplot2::aes_string(fill = "Group"), 
        stat = "identity"
      )
    ) +
    ggplot2::ylab("Inferred Ancestry Proportions") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
  if(label.pops) {
    g <- g + 
      ggplot2::geom_vline(xintercept = pop.div, size = 1.5) +
      ggplot2::scale_x_continuous(
        name = "", 
        breaks = pop.cntr, 
        labels = names(pop.cntr),
        expand = c(0, 0)
      )
  } else if (label.inds) {
    g <- g + 
      ggplot2::geom_vline(xintercept = pop.div, size = 1.5) +
      ggplot2::scale_x_continuous(
        name = "", 
        breaks = ind.cntr, 
        labels = names(ind.cntr),
        expand = c(0, 0)
      ) +
      ggplot2::theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  } else {
    g <- g + 
      ggplot2::xlab("") + 
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  if(horiz) g <- g + ggplot2::coord_flip()
  if(!is.null(col)) g <- g + ggplot2::scale_fill_manual(values = col)
  
  if(plot) print(g)
  invisible(g)
}



################################################################
###  Wrapper function around conStruct to do a whole bunch   ###
################################################################

conStructWrapper <- function(file = NULL, ind.version = F, pop.version = F, chains = NULL, iterations = NULL) {
  
  # Get model from file name
  model <- strsplit(basename(file), split = "\\.")[[1]][1]
  # Get generation from file name
  gen <- strsplit(basename(file), split = "\\.")[[1]][3]
  # Make prefix for output files
  prefix <- paste0(model, ".gen.", gen)
  print(prefix)
  
  # Info for helping format input files:
  #vignette(topic="format-data",package="conStruct")
  
  # Get allele frequencies for conStruct, from structure file
  allele.freqs <- structure2conStruct(infile = file,
                                      onerowperind = FALSE,
                                      start.loci = 2,
                                      start.samples = 2,
                                      missing.datum = -9,
                                      outfile = "temp")
  # Remove auto-generated output file
  unlink("temp.RData")
  
  
  
  ######### Individual version:
  
  if (ind.version) {
    
    # Make matrix of geographic sampling coordinates
    # x and y (made up) instead of lat and long
    geo.coords <- matrix(nrow = nrow(allele.freqs), ncol = 2)
    rownames(geo.coords) <- rownames(allele.freqs)
    colnames(geo.coords) <- c("x", "y")
    geo.coords[,"x"] <- sort(rep(0:20,20))
    geo.coords[,"y"] <- rep(1:20,21)
    
    # Make matrix of geographic distances
    geo.dists <- matrix(nrow = nrow(geo.coords), ncol = nrow(geo.coords))
    rownames(geo.dists) <- rownames(geo.coords)
    colnames(geo.dists) <- rownames(geo.dists)
    for (i in rownames(geo.dists)) {
      for (j in rownames(geo.dists)) {
        pop1 <- geo.coords[i,"x"]
        pop2 <- geo.coords[j,"x"]
        dist <- abs(pop2 - pop1)
        geo.dists[i,j] <- dist
        geo.dists[j,i] <- dist
      }
    }
  }
  
  
  if (pop.version) {
    ######### Population version:
    # Collapse allele frequencies into populations
    allele.freqs <- data.frame(cbind(sort(rep(0:20,20)), allele.freqs))
    pop.freqs <- matrix(nrow = 21, ncol = ncol(allele.freqs))
    for (i in 0:20) {
      pop.freqs[i+1,] <- colMeans(allele.freqs[allele.freqs$X1==i,])
    }
    # Get rid of pop column
    pop.freqs <- pop.freqs[,2:ncol(pop.freqs)]
    row.names(pop.freqs) <- 0:20
    
    
    # Make matrix of geographic sampling coordinates
    # x and y (made up) instead of lat and long
    geo.coords <- matrix(nrow = nrow(pop.freqs), ncol = 2)
    rownames(geo.coords) <- rownames(pop.freqs)
    colnames(geo.coords) <- c("x", "y")
    geo.coords[,"x"] <- 0:20
    geo.coords[,"y"] <- rep(1,21)
    
    # Make matrix of geographic distances
    geo.dists <- matrix(nrow = nrow(geo.coords), ncol = nrow(geo.coords))
    rownames(geo.dists) <- rownames(geo.coords)
    colnames(geo.dists) <- rownames(geo.dists)
    for (i in rownames(geo.dists)) {
      for (j in rownames(geo.dists)) {
        pop1 <- geo.coords[i,"x"]
        pop2 <- geo.coords[j,"x"]
        dist <- abs(pop2 - pop1)
        geo.dists[i,j] <- dist
        geo.dists[j,i] <- dist
        
      }
    }
    allele.freqs <- pop.freqs
  }
  
  
  # how to run a conStruct analysis
  #vignette(topic="run-conStruct",package="conStruct")
  
  # Run spatial conStruct for K=2
  run <- conStruct(spatial = TRUE, 
                   K = 2, 
                   freqs = allele.freqs, 
                   geoDist = geo.dists, 
                   coords = geo.coords, 
                   prefix = prefix, 
                   n.chains = chains, 
                   n.iter = iterations, 
                   make.figs = T, 
                   save.files = T)
  
  # Save the run to file
  save(run, file = paste0("conStructRuns/", prefix, ".run"))
  
  # how to visualize the output of a conStruct model
  #vignette(topic="visualize-results",package="conStruct")
  
  admix.props <- run$chain_1$MAP$admix.proportions
  #make.structure.plot(admix.proportions = admix.props, sample.names = 0:20)
  
  #make.admix.pie.plot(admix.proportions = admix.props,
  #                    coords = geo.coords,
  #                    radii = 3)
  
  color_ramp <- colorRampPalette(c("blue", "red"))
  cols <- color_ramp(20)
  
  
  qmat <- data.frame(cbind(0:20, rep(0,21), 0:20, admix.props))
  colnames(qmat) <- c("id", "pct.miss", "orig.pop", "Group.1", "Group.2")
  # make pop 0 blue
  if (qmat[1,"Group.1"] < 0.5) {
    qmat$Group.1 <- 1 - qmat$Group.1
    qmat$Group.2 <- 1 - qmat$Group.2
  }
  conStructPlot <- betterStructurePlot(qmat, label.pops = F, horiz = F, legend.position = "top", col = color_ramp(2), sort.probs = F, type = "bar")
  conStructPlot <- conStructPlot + labs(title = paste0(model, ", gen: ", gen))
  ggplot2::ggsave(filename = paste0("conStructPlots/", model, ".", gen, ".conStruct.png"), conStructPlot, width = 10, height = 5, path = getwd())
  
}

#######################################
####      End define functions     ####
#######################################


# Get arguments from command line
f = commandArgs(trailingOnly=TRUE)
print(f)

conStructWrapper(file = f, pop.version = T, chains = 5, iterations = 200000)



