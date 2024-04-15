setwd("/panfs/pfs.local/work/colella/ben/002.simulations/003.SecondaryContact_vs_IBD/03.conStruct/crossValidate")

library(ggplot2) 
library(conStruct)
library(scales)
library(reshape2)
library(doParallel)


################################################################
###  Wrapper function around conStruct to do a whole bunch   ###
################################################################

crossValWrapper <- function(file = NULL, ind.version = F, pop.version = F, iterations = NULL,
                            train.prop = NULL, n.reps = NULL, n.nodes = NULL, prefix = NULL) {
  
  # Get model from file name
  model <- strsplit(file, split = "\\.")[[1]][1]
  model <- strsplit(model, split = "/")[[1]][2]
  # Get generation from file name
  gen <- strsplit(file, split = "\\.")[[1]][3]
  
  
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
  
  # Run cross validation for spatial and nonspatial conStruct for K=1:3
  my.xvals <- x.validation(train.prop = train.prop,
                           n.reps = n.reps,
                           K = 1:3,
                           freqs = allele.freqs,
                           data.partitions = NULL,
                           geoDist = geo.dists,
                           coords = geo.coords,
                           prefix = prefix,
                           n.iter = iterations,
                           make.figs = FALSE,
                           save.files = FALSE,
                           parallel = TRUE,
                           n.nodes = n.nodes)
  
  
  # format results from the output list
  sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
  nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)
 
  
  # first, get the 95% confidence intervals for the spatial and nonspatial
  #   models over values of K (mean +/- 1.96 the standard error)
  
  sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
  nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
  
  # then, plot cross-validation results for K=1:3 with 8 replicates
  
  par(mfrow=c(1,2))
  plot(rowMeans(sp.results),
       pch=19,col="blue",
       ylab="predictive accuracy",xlab="values of K",
       ylim=range(sp.results,nsp.results),
       main="cross-validation results")
  points(rowMeans(nsp.results),col="green",pch=19)
  
  # finally, visualize results for the spatial model
  #   separately with its confidence interval bars
  #
  # note that you could do the same with the spatial model, 
  #   but the confidence intervals don't really show up 
  #   because the differences between predictive accuracies
  #   across values of K are so large.
  
  plot(rowMeans(sp.results),
       pch=19,col="blue",
       ylab="predictive accuracy",xlab="values of K",
       ylim=range(sp.CIs),
       main="spatial cross-validation results")
  segments(x0 = 1:nrow(sp.results),
           y0 = sp.CIs[1,],
           x1 = 1:nrow(sp.results),
           y1 = sp.CIs[2,],
           col = "blue",lwd=2)
  
   
}



# Get arguments from command line
ags = commandArgs(trailingOnly=TRUE)
f <- ags[1]
prefix <- ags[2]
print(f)

#### Run cross validation

crossValWrapper(file = f, pop.version = T, iterations = 50,
                  train.prop = 0.8, n.reps = 1, n.nodes = 8, prefix = prefix)

dev.print(pdf, paste0("prefix",".xVal.pdf"))



