setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")

library("PopGenReport")
library("adegenet")
library("tinytex")
library("poppr")
library("pegas")
library("diveRsity")
library("hierfstat")
library("mmod")
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap) 
library(diveRsity)
library(dartR)
library(RColorBrewer)
library(gtools)
library(reshape2)
library(ggplot2) 
library(cowplot)

library(dplyr)
library(apex)
library(pcadapt)
library(vecsets)
library(readr)
library(gridExtra)

library(SNPfiltR)
library(readxl)

library(strataG)
library(pophelper)
library(scales)

library(gifski)
library(gtools)

library(patchwork)

library(ggpubr)

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
    ggplot2::ylab("Anc. Pro.") +
    ggplot2::scale_y_continuous(n.breaks = 3, expand = c(0, 0)) +
    ggplot2::theme(
      axis.title.y = element_text(size = 30),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
  if(label.pops) {
    g <- g + 
      ggplot2::geom_vline(xintercept = pop.div, size = 1) +
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




############################################################
###    Define function to make PCAs and Triangle plots   ###
############################################################

publicationPlotter <- function(file = NULL, model = NULL, 
                              downsample.pops = NA, inds.per.pop = 20,
                              missing.perc = NA, missing.plot = F, missplots = F,
                              fixed.diffs.triangle = F, p1 = NULL, p2 = NULL, difference = NULL,
                              plot = T, pcatri = F, tri = F) {
  # Make R objects
  vcfR <- read.vcfR(file)
  gl <- vcfR2genlight(vcfR)
  gi <- vcfR2genind(vcfR)
  
  # Get generation from file name
  gen <- strsplit(file, split = "\\.")[[1]][2]
  
  # Make popmap and assign populations
  sim.popmap <- data.frame(matrix(nrow = length(colnames(extract.gt(vcfR)))))
  colnames(sim.popmap) <- "id"
  sim.popmap$id <- colnames(extract.gt(vcfR))
  num.pops <- nrow(sim.popmap) / 20
  sim.popmap$pop <- sort(c(rep(0:(num.pops-1),20)))
  pop(gl) <- sim.popmap$pop

  
  
  # Downsample pops if asked to
  if (!is.na(downsample.pops)) {
    sim.popmap <- sim.popmap[sim.popmap$pop %in% downsample.pops,]
    vcfR <- vcfR[samples = sim.popmap$id]
    gl <- vcfR2genlight(vcfR)
    pop(gl) <- sim.popmap$pop
  }
  
  
  # Downsample inds if asked to:
  if (!inds.per.pop==20) {
    temp.popmap <- data.frame(matrix(nrow=0,ncol=ncol(sim.popmap)))
    for (a in unique(sim.popmap$pop)) {
      temp.popmap <- rbind(temp.popmap, sim.popmap[sim.popmap$pop==a,][1:inds.per.pop,])
    }
    sim.popmap <- temp.popmap
    vcfR <- vcfR[samples = sim.popmap$id]
    gl <- vcfR2genlight(vcfR)
    pop(gl) <- sim.popmap$pop
  }
  
  
  # Add missing data if asked to 
  if (!is.na(missing.perc)) {
    # use beta distribution, alpha always 1.5, beta is whatever gives the overall percentage of missing data
    b <- 1.5 * ((1/missing.perc) - 1)
    miss <- rbeta(ncol(extract.gt(vcfR)), shape1 = 1.5, shape2 = b)
    nsnps <- nrow(vcfR@gt)
    miss <- round(miss * nsnps)
    names(miss) <- colnames(vcfR@gt)[2:length(colnames(vcfR@gt))]
    for (i in names(miss)) {
      indices <- sample(1:nsnps, size = miss[i], replace = F)
      vcfR@gt[indices,i] <- NA
    }
    gl <- vcfR2genlight(vcfR)
    pop(gl) <- sim.popmap$pop
    gi <- vcfR2genind(vcfR)
  }
  
  
  if (plot) {
    #perform PCA
    pca<-glPca(gl, nf = 10)
    #pull pca scores out of df
    pca.scores<-as.data.frame(pca$scores)
    pca.scores$id <- sim.popmap$id
    pca.scores$pop <- sim.popmap$pop
    #calculate missingness by individual
    ms<-colSums(is.na(vcfR@gt))/nrow(vcfR@gt)
    #add missingness to df
    pca.scores$missing<-ms[-1]
    #record variance percentages explained
    var.frac<- pca$eig/sum(pca$eig)*100
    # Flip PC1 if necessary
    ave.subpop.p1 <- mean(pca.scores[pca.scores$pop==p1,"PC1"])
    ave.subpop.p2 <- mean(pca.scores[pca.scores$pop==p2,"PC1"])
    if (ave.subpop.p1 > ave.subpop.p2) {
      pca.scores$PC1 <- pca.scores$PC1 * -1
    }    
    # Flip PC2 if necessary
    ave.subpop.p1 <- mean(pca.scores[pca.scores$pop==p1,"PC2"])
    median.pop <- which.min(abs(unique(pca.scores$pop) - floor(median(unique(pca.scores$pop)))))
    median.pop <- unique(pca.scores$pop)[median.pop]
    ave.subpop.median <- mean(pca.scores[pca.scores$pop==median.pop,"PC2"])
    if (ave.subpop.p1 > ave.subpop.median) {
      pca.scores$PC2 <- pca.scores$PC2 * -1
    }

    
    
    # put dataframe in different order for plotting purposes
    pca.scores <- rbind(pca.scores[pca.scores$pop %in% 10:20,], pca.scores[pca.scores$pop %in% 0:9,][order(nrow(pca.scores[pca.scores$pop %in% 0:9,]):1),] )
    
    
    #ggplot color by pop
    pca.plot <- ggplot(pca.scores, aes(x=PC1, y=PC2, color=pop)) +
      geom_point(cex = 2, alpha=1)+
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab("PC1")+
      ylab("PC2")+
      #xlab(paste("PC1, ", sprintf('%.2f', var.frac[1]), "% var. exp.", sep = ""))+
      #ylab(paste("PC2, ", sprintf('%.2f', var.frac[2]), "% var. exp.", sep = ""))+
      #labs(title = paste0(model, ", generation: ", gen, ", total missing data: ", missing.perc)) +
      scale_color_gradient2(low = "#313695", mid = "yellow3", high = "#a50026", midpoint = 10) +
      theme_classic() +
      theme(legend.position = "none", axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 30),
            axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))

    #print(pca.plot)

    
    
    if (missing.plot) {
      #ggplot color by missing
      pca.missing <- ggplot(pca.scores, aes(x=PC1, y=PC2, color=missing)) +
        geom_point(cex = 2, alpha=1)+
        guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
        xlab("PC1")+
        ylab("PC2")+
        #xlab(paste("PC1, ", sprintf('%.2f', var.frac[1]), "% var. exp.", sep = ""))+
        #ylab(paste("PC2, ", sprintf('%.2f', var.frac[2]), "% var. exp.", sep = ""))+
        #labs(title = paste0(model, ", generation: ", gen, ", total missing data: ", missing.perc)) +
        scale_color_gradient(low = "#000000", high = "#f0f0f0", limits = c(0,1)) +
        theme_classic() +
        theme(legend.position = "none", axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 30),
              axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      #print(pca.missing)
      
      
    }
    

    
    
    if (fixed.diffs.triangle) {
      
      # get differences above allele frequency threshold
      vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = p1, p2 = p2, difference = difference)
      # get hybrid index
      hi <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = p1, p2 = p2)
      
      # number of differences above threshold
      d <- nrow(vcfR.diff@gt)
      # number of total SNPs
      t <- nrow(vcfR@gt)
      # percent of SNPs of threshold
      perc_snps <- round((d/t)*100, 1)
      
      # add id column
      hi$id <- rownames(hi)

      # get heterozygosity (don't count NAs as FALSE, leave as NA)
      het <- is.het(extract.gt(vcfR.diff), na_is_false = F)
      for (ind in sim.popmap$id) {
        # divide the number of sites that are het by the number of sites that have data
        hi[hi$id == ind, "het"] <- sum(het[,ind], na.rm = T)/sum(!is.na(het[,ind]))
      }
      # add pop to hybrid index dataframe
      for (ind in sim.popmap$id) {
        hi[hi$id == ind, "pop"] <- sim.popmap[sim.popmap$id == ind, "pop"]
      }
      # add missingness to hybrid index dataframe
      for (ind in pca.scores$id) {
        hi[hi$id == ind, "missing"] <- pca.scores[pca.scores$id == ind, "missing"]
      }
      
      # put dataframe in different order for plotting purposes
      hi <- rbind(hi[hi$pop %in% 10:20,], hi[hi$pop %in% 0:9,][order(nrow(hi[hi$pop %in% 0:9,]):1),] )
      
      #ggplot color by pop
      diff.triangle <- ggplot(hi, aes(x=hybrid.index, y=het, color=pop)) +
        stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
        geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
        geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
        geom_jitter(cex = 2, alpha=1, width = 0.02, height = 0.02)+
        guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
        xlab(paste("Hybrid Index"))+
        ylab(paste("Heterozygosity"))+
        #labs(title = paste0(model, ", ", d, " Differences (threshold: ", difference, "), generation: ", gen, ", total missing data: ", missing.perc)) +
        scale_color_gradient2(low = "#313695", mid = "yellow3", high = "#a50026", midpoint = 10) +
        #ylim(c(-0.05,1.05)) +
        #xlim(c(-0.05,1.05)) +
        scale_x_continuous(n.breaks = 3) + 
        scale_y_continuous(n.breaks = 3) +
        theme_classic() +
        theme(legend.position = "none", axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 30),
              axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      #theme(legend.position = c(0.13,0.75))
      
      #print(diff.triangle)

      diff.triangle.miss <- ggplot(hi, aes(x=hybrid.index, y=het, color=missing)) +
        stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
        geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
        geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
        geom_jitter(cex = 2, alpha=1, width = 0.02, height = 0.02)+
        guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
        xlab(paste("Hybrid Index"))+
        ylab(paste("Heterozygosity"))+
        #labs(title = paste0(model, ", ", d, " Differences (threshold: ", difference, "), generation: ", gen, ", total missing data: ", missing.perc)) +
        scale_color_gradient(low = "#000000", high = "#f0f0f0", limits = c(0,1)) +
        #ylim(c(-0.05,1.05)) +
        #xlim(c(-0.05,1.05)) +
        scale_x_continuous(n.breaks = 3) + 
        scale_y_continuous(n.breaks = 3) +
        theme_classic() +
        theme(legend.position = "none", axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 30),
              axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      
      #print(diff.triangle.miss)
    }
    
  }
  if (pcatri) {
    return(pca.plot + diff.triangle)
  }
  if (missplots) {
    return(pca.missing + diff.triangle.miss)
  }
  if (tri) {
    print(diff.triangle + annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0(perc_snps,"% of SNPs"), size=4))
    return(diff.triangle + annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0(perc_snps,"% of SNPs"), size=4))
  }

}




alleleFreqDiff <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  # make matrices of genotypes for each parental pop
  p1.gts <- data.frame(extract.gt(vcfR[samples = pm[pm$pop == p1,"id"]]))
  p2.gts <- data.frame(extract.gt(vcfR[samples = pm[pm$pop == p2,"id"]]))
  
  # make dataframe to keep track of frequency of "1" allele at each locus
  af <- data.frame(matrix(nrow = nrow(p1.gts), ncol = 2))
  colnames(af) <- c("p1", "p2")
  rownames(af) <- rownames(p1.gts)
  
  # recode missing data
  p1.gts[is.na(p1.gts)] <- -9
  p2.gts[is.na(p2.gts)] <- -9
  # recode alelles
  p1.gts[p1.gts=="0|0"] <- 0
  p1.gts[p1.gts=="0|1"] <- 1
  p1.gts[p1.gts=="1|0"] <- 1
  p1.gts[p1.gts=="1|1"] <- 2
  p1.gts[p1.gts=="0/0"] <- 0
  p1.gts[p1.gts=="0/1"] <- 1
  p1.gts[p1.gts=="1/0"] <- 1
  p1.gts[p1.gts=="1/1"] <- 2
  p2.gts[p2.gts=="0|0"] <- 0
  p2.gts[p2.gts=="0|1"] <- 1
  p2.gts[p2.gts=="1|0"] <- 1
  p2.gts[p2.gts=="1|1"] <- 2
  p2.gts[p2.gts=="0/0"] <- 0
  p2.gts[p2.gts=="0/1"] <- 1
  p2.gts[p2.gts=="1/0"] <- 1
  p2.gts[p2.gts=="1/1"] <- 2

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  for (i in 1:nrow(af)) {
    # p1:
    # get genotypes for locus, remove missing data
    loc.p1 <- p1.gts[i,]
    loc.p1 <- loc.p1[loc.p1!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p1"] <- sum(loc.p1)/(length(loc.p1)*2)
    
    #p2:
    # get genotypes for locus, remove missing data
    loc.p2 <- p2.gts[i,]
    loc.p2 <- loc.p2[loc.p2!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p2"] <- sum(loc.p2)/(length(loc.p2)*2)
  }
  # add allele frequency differences to af
  af$diff <- abs(af[,"p1"] - af[,"p2"])
  
  # get names of loci with allele frequency difference above threshold
  loci <- rownames(af[af$diff >= difference,])
  
  # get names of all loci in vcfR object
  all.loci <- rownames(extract.gt(vcfR))
  
  # get indices in vcfR object of loci with difference above theshold
  loci.indices <- all.loci %in% loci
  loci.indices <- which(loci.indices)
  
  # subset vcfR by fixed diff loci
  vcfR.diff <- vcfR[loci.indices]
  return(vcfR.diff)
}



hybridIndex <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL) {

  m <- extract.gt(vcfR)
  # recode missing data
  m[is.na(m)] <- -9
  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2
  # make new matrix of same size as m
  n <- matrix(nrow = nrow(m), ncol = ncol(m))
  
  # make matrices of genotypes for each parental pop
  p1.gts <- data.frame(m[,pm[pm$pop == p1,"id"]])
  p2.gts <- data.frame(m[,pm[pm$pop == p2,"id"]])
  
  # make dataframe to keep track of frequency of "1" allele at each locus
  af <- data.frame(matrix(nrow = nrow(m), ncol = 2))
  colnames(af) <- c("p1", "p2")
  rownames(af) <- rownames(m)
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  for (i in 1:nrow(m)) {

    # p1:
    # get genotypes for locus, remove missing data
    loc.p1 <- p1.gts[i,]
    loc.p1 <- loc.p1[loc.p1!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p1"] <- sum(loc.p1)/(length(loc.p1)*2)
    
    #p2:
    # get genotypes for locus, remove missing data
    loc.p2 <- p2.gts[i,]
    loc.p2 <- loc.p2[loc.p2!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p2"] <- sum(loc.p2)/(length(loc.p2)*2)
    
    # for each locus, if parental pop 1 has a higher frequency of "2" allele, assign p1.allele as 2
    if (af[i,"p1"] > af[i,"p2"]) {
      p1.allele <- 2
      # else, if parental pop1 has a lower frequency of "2 allele", assign p1.allele as 0
    } else {
      p1.allele <- 0
    }
    
    # compare every individual to parental pop 1, giving scores of:
    # 0 = matching p1
    # 1 = being a het
    # 2 = not matching p1 (matching p2)
    # -9 = NA
    for (j in 1:ncol(m)) {
      if (m[i,j] == -9) {
        next
      }
      if (m[i,j] == p1.allele) {
        n[i,j] <- 0
      } else if (m[i,j] == 1) {
        n[i,j] <- 1
      } else {
        n[i,j] <- 2
      }
    }
  }
  colnames(n) <- colnames(m)
  rownames(n) <- rownames(m)
  # count alleles, removing NAs
  counts <- colSums(n, na.rm = T)
  # count nonmissing genotypes for each ind
  sites <- colSums(!apply(n, MARGIN = 2, is.na))
  # calculate hybrid index
  hi <- counts / (sites*2)
  hi <- data.frame(hi)
  colnames(hi) <- "hybrid.index"
  return(hi)
 
}





################################################################
#    Define function to make and save structure plot images    #
################################################################
StructureImage <- function(dir = NULL, files = NULL, model.name = NULL, gen = NULL, miss = NA) {
  
  # Import STRUCTURE results
  for (file in files) {
    # Get number to name structure run
    num_of_k <- paste(sapply(strsplit(file, split='_', fixed=TRUE), `[`, 2), sapply(strsplit(file, split='_', fixed=TRUE), `[`, 3), sep = "_")
    # read in structure results
    raw <- structureRead(paste0(dir,"/",file))
    # Make popmap
    sim.popmap <- data.frame(matrix(nrow = nrow(raw$q.mat)))
    colnames(sim.popmap) <- "id"
    sim.popmap$id <- raw$q.mat$id
    num.pops <- nrow(sim.popmap) / 20
    sim.popmap$pop <- sort(c(rep(0:(num.pops-1),20)))
    # Assign pop names
    for (i in 1:nrow(raw$q.mat)) {
      current_ind <- raw$q.mat[i,"id"]
      raw$q.mat[i,"orig.pop"] <- sim.popmap[sim.popmap$id == current_ind,"pop"]
    }
    # Rename to structure run
    assign(paste0(num_of_k, "_structure"), raw)
  }
  
  # if K2_rep1_structure is weird (pop 0 and pop 20 the same), use K2_rep2_structure
  if (K2_rep1_structure$q.mat[1,"Group.1"] > 0.95 & K2_rep1_structure$q.mat[420,"Group.1"] > 0.95) {
    K2_rep1_structure <- K2_rep2_structure
  }
  
  
  # make pop 0 blue
  if (K2_rep1_structure$q.mat[1,"Group.1"] < 0.5) {
    K2_rep1_structure$q.mat$Group.1 <- 1 - K2_rep1_structure$q.mat$Group.1
    K2_rep1_structure$q.mat$Group.2 <- 1 - K2_rep1_structure$q.mat$Group.2
  }
  
  #for (k in ls(pattern = "K*_structure")) {
  for (k in ls(pattern = "K2_rep1_structure")) {
    # assign temp object for data
    kdata <- get(k)
    num <- sapply(strsplit(paste(sapply(strsplit(k, split='_', fixed=TRUE), `[`, 1), sep = "_"), split = "K", fixed =TRUE), `[`, 2)
    plot <- betterStructurePlot(kdata$q.mat, label.pops = T, horiz = F, legend.position = "none", col = c("#313695", "#a50026"), sort.probs = T, type = "bar")
    #if (!is.na(miss)) {
    #  plot <- plot + labs(title = paste0(model.name, ", gen: ", gen, ", miss: ", miss))
    #} else {
    #  plot <- plot + labs(title = paste0(model.name, ", gen: ", gen))
    #}
  }
  return(plot)
}




############################################################
###    Define function to make PCAs and Triangle plots   ###
############################################################

plotsForGifs <- function(file = NULL, model = NULL, 
                               downsample.pops = NA, inds.per.pop = 20,
                               missing.perc = NA, missing.plot = F, missplots = F,
                               p1 = NULL, p2 = NULL, difference = NULL,
                               tri = F) {
  # Make R objects
  vcfR <- read.vcfR(file)
  
  # Get generation from file name
  gen <- strsplit(file, split = "\\.")[[1]][2]
  
  # Make popmap and assign populations
  sim.popmap <- data.frame(matrix(nrow = length(colnames(extract.gt(vcfR)))))
  colnames(sim.popmap) <- "id"
  sim.popmap$id <- colnames(extract.gt(vcfR))
  num.pops <- nrow(sim.popmap) / 20
  sim.popmap$pop <- sort(c(rep(0:(num.pops-1),20)))
  
  
  
  # Downsample pops if asked to
  if (!is.na(downsample.pops)) {
    sim.popmap <- sim.popmap[sim.popmap$pop %in% downsample.pops,]
    vcfR <- vcfR[samples = sim.popmap$id]
  }
  
  
  # Downsample inds if asked to:
  if (!inds.per.pop==20) {
    temp.popmap <- data.frame(matrix(nrow=0,ncol=ncol(sim.popmap)))
    for (a in unique(sim.popmap$pop)) {
      temp.popmap <- rbind(temp.popmap, sim.popmap[sim.popmap$pop==a,][1:inds.per.pop,])
    }
    sim.popmap <- temp.popmap
    vcfR <- vcfR[samples = sim.popmap$id]
  }
  
  
  # Add missing data if asked to 
  if (!is.na(missing.perc)) {
    # use beta distribution, alpha always 1.5, beta is whatever gives the overall percentage of missing data
    b <- 1.5 * ((1/missing.perc) - 1)
    miss <- rbeta(ncol(extract.gt(vcfR)), shape1 = 1.5, shape2 = b)
    nsnps <- nrow(vcfR@gt)
    miss <- round(miss * nsnps)
    names(miss) <- colnames(vcfR@gt)[2:length(colnames(vcfR@gt))]
    for (i in names(miss)) {
      indices <- sample(1:nsnps, size = miss[i], replace = F)
      vcfR@gt[indices,i] <- NA
    }
  }
  
  
  
  # get differences above allele frequency threshold
  vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = p1, p2 = p2, difference = difference)
  # get hybrid index
  hi <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = p1, p2 = p2)
  
  # number of differences above threshold
  d <- nrow(vcfR.diff@gt)
  # number of total SNPs
  t <- nrow(vcfR@gt)
  # percent of SNPs of threshold
  perc_snps <- round((d/t)*100, 1)
  
  # add id column
  hi$id <- rownames(hi)
  
  # get heterozygosity (don't count NAs as FALSE, leave as NA)
  het <- is.het(extract.gt(vcfR.diff), na_is_false = F)
  for (ind in sim.popmap$id) {
    # divide the number of sites that are het by the number of sites that have data
    hi[hi$id == ind, "het"] <- sum(het[,ind], na.rm = T)/sum(!is.na(het[,ind]))
  }
  # add pop to hybrid index dataframe
  for (ind in sim.popmap$id) {
    hi[hi$id == ind, "pop"] <- sim.popmap[sim.popmap$id == ind, "pop"]
  }
  # add missingness to hybrid index dataframe
  #calculate missingness by individual
  ms<-colSums(is.na(vcfR@gt))/nrow(vcfR@gt)
  ms<-ms[-1]
  
  for (ind in names(ms)) {
    hi[hi$id == ind, "missing"] <- ms[ind]
  }
  
  # put dataframe in different order for plotting purposes
  hi <- rbind(hi[hi$pop %in% 10:20,], hi[hi$pop %in% 0:9,][order(nrow(hi[hi$pop %in% 0:9,]):1),] )
  
  #ggplot color by pop
  diff.triangle <- ggplot(hi, aes(x=hybrid.index, y=het, color=pop)) +
    stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
    geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
    geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
    geom_jitter(cex = 2, alpha=1, width = 0.02, height = 0.02)+
    guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
    xlab(paste("Hybrid Index"))+
    ylab(paste("Heterozygosity"))+
    #labs(title = paste0(model, ", ", d, " Differences (threshold: ", difference, "), generation: ", gen, ", total missing data: ", missing.perc)) +
    scale_color_gradient2(low = "#313695", mid = "yellow3", high = "#a50026", midpoint = 10) +
    ylim(c(-0.05,1.05)) +
    xlim(c(-0.05,1.05)) +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))
  #theme(legend.position = c(0.13,0.75))
  
  #print(diff.triangle)
  
  diff.triangle.miss <- ggplot(hi, aes(x=hybrid.index, y=het, color=missing)) +
    stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
    geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
    geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
    geom_jitter(cex = 2, alpha=1, width = 0.02, height = 0.02)+
    guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
    xlab(paste("Hybrid Index"))+
    ylab(paste("Heterozygosity"))+
    #labs(title = paste0(model, ", ", d, " Differences (threshold: ", difference, "), generation: ", gen, ", total missing data: ", missing.perc)) +
    scale_color_gradient(low = "#000000", high = "#f0f0f0", limits = c(0,1)) +
    ylim(c(-0.05,1.05)) +
    xlim(c(-0.05,1.05)) +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))
  
  #print(diff.triangle.miss)

  if (missplots) {
    return(diff.triangle.miss)
  }
  if (tri) {
    print(diff.triangle + annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0("gen: ", gen), size=6))
    return(diff.triangle + annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0("gen: ", gen), size=6))
  }
  
}



######################################################################
####     End define functions   #####
######################################################################


#  Fig 1 triangle plot schematic
example.triangle <- data.frame(cbind(c(0,0.25,0.5,0.5,0.75,1,0.5/4,0.5/8,0.5/16,1-0.5/4,1-0.5/8,1-0.5/16,-0.1,1.1,1),c(0,0.5,1,0.5,0.5,0,1/4,1/8,1/16,1/4,1/8,1/16,0.5,0.5,1.1)))
colnames(example.triangle) <- c("hi", "het")
example.triangle$class <- c("P1", "BC1", "F1", "F2", "BC2", "P2", "BC1.2", "BC1.3", "BC1.4", "BC2.2", "BC2.3", "BC2.4", "z", "z", "z")


example.triangle.plot <- ggplot(example.triangle, aes(x=hi, y=het, color = class)) +
  stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
  geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
  geom_point(cex = 3) +
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("Hybrid Index"))+
  ylab(paste("Heterozygosity"))+
  scale_color_manual(values = c("#e7d4e8", "#c2a5cf", "#9970ab", "#762a83",
                                "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837",  
                                "#bababa", "#bababa", "#40004b", "#00441b",
                                "#FFFFFF", "#FFFFFF", "#FFFFFF")) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_classic() +
  annotate("text", x=-0.075,y=-0.02,label="P1",color="#000000") +
  annotate("text", x=1.075,y=-0.02,label="P2",color="#000000") +
  annotate("text", x=0.5,y=1.075,label="F1",color="#000000") +
  annotate("text", x=0.5,y=0.575,label="F2 and later",color="#000000") +
  annotate("text", x=0.155,y=0.5,label="BC1",color="#000000") +
  annotate("text", x=0.845,y=0.5,label="BC1",color="#000000") +
  annotate("text", x=0.03,y=0.25,label="BC2",color="#000000") +
  annotate("text", x=-0.032,y=0.155,label="BC3",color="#000000") +
  annotate("text", x=-0.06075,y=0.0725,label="BC4",color="#000000") +
  annotate("text", x=0.97,y=0.25,label="BC2",color="#000000") +
  annotate("text", x=1.032,y=0.155,label="BC3",color="#000000") +
  annotate("text", x=1.06075,y=0.0725,label="BC4",color="#000000") +
  theme(legend.position = "none")


setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Figure1.pdf", example.triangle.plot, 
                width = 3.75, height = 2.5, path = getwd(), limitsize = F)

#   Fig 2 is schematic of simulations


##########################
#        Figure 3        #
##########################

# PCA and triangle plots
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")
for (model in c("IBD", "SC", "TZ")) {
  if (model == "IBD") {
    full.model <- "IBD"
  }
  if (model == "SC") {
    full.model <- "SecondaryContact"
  }
  if (model == "TZ") {
    full.model <- "TensionZone"
  }
  for (f in c(0,200,1000,10000,19000)) {
    patch_design <- c(area(1, 1, 1, 4),          # Specify design of grid
                      area(2, 1, 3, 4))
    if (exists("fig3")) {
      temp <- publicationPlotter(file = paste0(full.model,"/gen.",f,".vcf"), pcatri = T, fixed.diffs.triangle = T, difference = 0.5, p1 = 0, p2 = 20)
      dir <- paste0(full.model,"/str.output/01.output.",model,".gen.",f,".200K.K2/")
      files <- list.files(dir)
      files <- files[grepl("*f$", files)]  
      temp2 <- StructureImage(dir = dir, files = files, model.name = model, gen = f)
      temp2 <- temp2 + theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      temp <- temp2 / temp + plot_layout(design = patch_design)
      fig3 <- fig3 / temp
    } else {
      fig3 <- publicationPlotter(file = paste0(full.model,"/gen.",f,".vcf"), pcatri = T, fixed.diffs.triangle = T, difference = 0.5, p1 = 0, p2 = 20)
      dir <- paste0(full.model,"/str.output/01.output.",model,".gen.",f,".200K.K2/")
      files <- list.files(dir)
      files <- files[grepl("*f$", files)]  
      temp <- StructureImage(dir = dir, files = files, model.name = model, gen = f)
      temp <- temp + theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      fig3 <- temp / fig3 + plot_layout(design = patch_design)
    }
  }

  assign(paste0(model, ".fig3"), fig3)
  rm(fig3)
}


fig3 <- ggarrange(IBD.fig3[[1]]/IBD.fig3[[2]]+ plot_layout(design = patch_design)
                  , SC.fig3[[1]]/SC.fig3[[2]]+ plot_layout(design = patch_design)
                  , TZ.fig3[[1]]/TZ.fig3[[2]]+ plot_layout(design = patch_design)
                  , 
          IBD.fig3[[3]], SC.fig3[[3]], TZ.fig3[[3]], 
          IBD.fig3[[4]], SC.fig3[[4]], TZ.fig3[[4]], 
          IBD.fig3[[5]], SC.fig3[[5]], TZ.fig3[[5]], 
          IBD.fig3[[6]], SC.fig3[[6]], TZ.fig3[[6]], 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
          font.label = list(size = 30),
          #label.y = 0.1,
          #label.x = 0.1,
          ncol = 3, nrow = 5)



setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Figure3_fromR_test.pdf", fig3, 
                width = 30, height = 33, path = getwd(), limitsize = F)


fig3.small <- ggarrange(IBD.fig3[[1]]/IBD.fig3[[2]]+ plot_layout(design = patch_design)
                  , SC.fig3[[1]]/SC.fig3[[2]]+ plot_layout(design = patch_design)
                  , TZ.fig3[[1]]/TZ.fig3[[2]]+ plot_layout(design = patch_design)
                  , 
                  IBD.fig3[[4]], SC.fig3[[4]], TZ.fig3[[4]], 
                  IBD.fig3[[5]], SC.fig3[[5]], TZ.fig3[[5]], 
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                  font.label = list(size = 30),
                  #label.y = 0.1,
                  #label.x = 0.1,
                  ncol = 3, nrow = 3)

ggplot2::ggsave(filename = "Figure3_fromR_small.pdf", fig3.small, 
                width = 30, height = 20, path = getwd(), limitsize = F)



##########################
#        Figure 4        #
##########################

# PCA and triangle plots
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")
for (model in c("IBD", "SC", "TZ")) {
  if (model == "IBD") {
    full.model <- "IBD"
  }
  if (model == "SC") {
    full.model <- "SecondaryContact"
  }
  if (model == "TZ") {
    full.model <- "TensionZone"
  }
  for (m in c(0.1,0.2,0.3,0.4,0.5)) {
    patch_design <- c(area(1, 1, 1, 4),          # Specify design of grid
                      area(2, 1, 3, 4))
    if (exists("fig4")) {
      temp <- publicationPlotter(file = paste0(full.model,"/gen.1000.vcf"), fixed.diffs.triangle = T, difference = 0.5, p1 = 0, p2 = 20, missing.perc = m, missing.plot = T, missplots = T)
      dir <- paste0(full.model,"/str.miss.output/01.output.",model,".gen.1000.miss.",m,".200K.K2/")
      files <- list.files(dir)
      files <- files[grepl("*f$", files)]  
      temp2 <- StructureImage(dir = dir, files = files, model.name = model, gen = "1000")
      temp2 <- temp2 + theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      temp <- temp2 / temp + plot_layout(design = patch_design)
      fig4 <- fig4 / temp
    } else {
      fig4 <- publicationPlotter(file = paste0(full.model,"/gen.1000.vcf"), fixed.diffs.triangle = T, difference = 0.5, p1 = 0, p2 = 20, missing.perc = m, missing.plot = T, missplots = T)
      dir <- paste0(full.model,"/str.miss.output/01.output.",model,".gen.1000.miss.",m,".200K.K2/")
      files <- list.files(dir)
      files <- files[grepl("*f$", files)]  
      temp <- StructureImage(dir = dir, files = files, model.name = model, gen = "1000") 
      temp <- temp + theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))
      fig4 <- temp / fig4 + plot_layout(design = patch_design)
      
    }
  }

  assign(paste0(model, ".fig4"), fig4)
  rm(fig4)
}




fig4 <- ggarrange(IBD.fig4[[1]]/IBD.fig4[[2]]+ plot_layout(design = patch_design)
                  , SC.fig4[[1]]/SC.fig4[[2]]+ plot_layout(design = patch_design)
                  , TZ.fig4[[1]]/TZ.fig4[[2]]+ plot_layout(design = patch_design)
                  , 
                  IBD.fig4[[3]], SC.fig4[[3]], TZ.fig4[[3]], 
                  IBD.fig4[[4]], SC.fig4[[4]], TZ.fig4[[4]], 
                  IBD.fig4[[5]], SC.fig4[[5]], TZ.fig4[[5]], 
                  IBD.fig4[[6]], SC.fig4[[6]], TZ.fig4[[6]], 
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
                  font.label = list(size = 30),
                  #label.y = 0.1,
                  #label.x = 0.1,
                  ncol = 3, nrow = 5)



setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Figure4_fromR_test.pdf", fig4, 
                width = 30, height = 33, path = getwd(), limitsize = F)


fig4.small <- ggarrange(IBD.fig4[[1]]/IBD.fig4[[2]]+ plot_layout(design = patch_design)
                  , SC.fig4[[1]]/SC.fig4[[2]]+ plot_layout(design = patch_design)
                  , TZ.fig4[[1]]/TZ.fig4[[2]]+ plot_layout(design = patch_design)
                  , 
                  IBD.fig4[[4]], SC.fig4[[4]], TZ.fig4[[4]], 
                  IBD.fig4[[6]], SC.fig4[[6]], TZ.fig4[[6]], 
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                  font.label = list(size = 30),
                  #label.y = 0.1,
                  #label.x = 0.1,
                  ncol = 3, nrow = 3)

ggplot2::ggsave(filename = "Figure4_fromR_small.pdf", fig4.small, 
                width = 30, height = 20, path = getwd(), limitsize = F)


##########################
#        Figure 5        #
##########################

# PCA and triangle plots
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")
for (model in c("IBD", "SC", "TZ")) {
  if (model == "IBD") {
    full.model <- "IBD"
  }
  if (model == "SC") {
    full.model <- "SecondaryContact"
  }
  if (model == "TZ") {
    full.model <- "TensionZone"
  }
  for (g in c(1000,10000)) {
    for (a in c(0.5,0.75,1)) {
      
      if (exists("fig5")) {
        temp <- publicationPlotter(file = paste0(full.model,"/gen.",g,".vcf"), fixed.diffs.triangle = T, difference = a, p1 = 0, p2 = 20, tri = T)
        fig5 <- fig5 / temp
      } else {
        fig5 <- publicationPlotter(file = paste0(full.model,"/gen.",g,".vcf"), fixed.diffs.triangle = T, difference = a, p1 = 0, p2 = 20, tri = T)
        
      }
    }
  }


  assign(paste0(model, ".fig5"), fig5)
  rm(fig5)
}


fig5 <- ggarrange(IBD.fig5[[1]], SC.fig5[[1]], TZ.fig5[[1]], 
          IBD.fig5[[2]], SC.fig5[[2]], TZ.fig5[[2]], 
          IBD.fig5[[3]], SC.fig5[[3]], TZ.fig5[[3]], 
          IBD.fig5[[4]], SC.fig5[[4]], TZ.fig5[[4]], 
          IBD.fig5[[5]], SC.fig5[[5]], TZ.fig5[[5]], 
          IBD.fig5[[6]], SC.fig5[[6]], TZ.fig5[[6]], 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R"),
          font.label = list(size = 14),
          #label.y = -0.2,
          label.x = 0.1,
          ncol = 3, nrow = 6)



setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Figure5_fromR.pdf", fig5, 
                width = 9, height = 15, path = getwd(), limitsize = F)





##########################
#        Figure 6        #
##########################

# PCA and triangle plots
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")
for (model in c("IBD", "SC", "TZ")) {
  if (model == "IBD") {
    full.model <- "IBD"
  }
  if (model == "SC") {
    full.model <- "SecondaryContact"
  }
  if (model == "TZ") {
    full.model <- "TensionZone"
  }
  for (n in c(20,10,5,2)) {
    
    if (exists("fig6")) {
      temp <- publicationPlotter(file = paste0(full.model,"/gen.1000.vcf"), fixed.diffs.triangle = T, difference = 0.5, p1 = 0, p2 = 20, tri = T, inds.per.pop = n)
      fig6 <- fig6 / temp
    } else {
      fig6 <- publicationPlotter(file = paste0(full.model,"/gen.1000.vcf"), fixed.diffs.triangle = T, difference = 0.5, p1 = 0, p2 = 20, tri = T, inds.per.pop = n)
      
    }
  }

  assign(paste0(model, ".fig6"), fig6)
  rm(fig6)
}



fig6 <- ggarrange(IBD.fig6[[1]], SC.fig6[[1]], TZ.fig6[[1]], 
                  IBD.fig6[[2]], SC.fig6[[2]], TZ.fig6[[2]], 
                  IBD.fig6[[3]], SC.fig6[[3]], TZ.fig6[[3]], 
                  IBD.fig6[[4]], SC.fig6[[4]], TZ.fig6[[4]], 
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                  font.label = list(size = 14),
                  #label.y = -0.2,
                  label.x = 0.1,
                  ncol = 3, nrow = 4)



setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Figure6_fromR.pdf", fig6, 
                width = 9, height = 10, path = getwd(), limitsize = F)


# Make color scales
dummy <- data.frame(x=c(0:20),pop=c(0:20),percent.missing=c(0:20))
byr.scale <- ggplot(dummy, aes(x=x, y=pop, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  scale_color_gradient2(low = "#313695", mid = "yellow3", high = "#a50026", midpoint = 10, n.breaks = 21) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(6, 'cm')) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size=30))

bw.scale <- ggplot(dummy, aes(x=x, y=pop, color=percent.missing)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  scale_color_gradient(low = "#000000", high = "#f0f0f0", limits = c(0,1)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(6, 'cm')) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size=30))

scale <- ggarrange(byr.scale, bw.scale,
                   ncol = 1)

setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "scales.from.R.pdf", scale, 
                width = 30, height = 10, path = getwd(), limitsize = F)

byr.scale.small <- ggplot(dummy, aes(x=x, y=pop, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  scale_color_gradient2(low = "#313695", mid = "yellow3", high = "#a50026", midpoint = 10, n.breaks = 21) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(3, 'cm')) +
  theme(legend.title = element_text(size=14), legend.text = element_text(size=14))

bw.scale.small <- ggplot(dummy, aes(x=x, y=pop, color=percent.missing)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  scale_color_gradient(low = "#000000", high = "#f0f0f0", limits = c(0,1)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(3, 'cm')) +
  theme(legend.title = element_text(size=14), legend.text = element_text(size=14))

scale.small <- ggarrange(byr.scale.small, bw.scale.small,
                   ncol = 1)

setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "small.scales.from.R.pdf", scale.small, 
                width = 30, height = 10, path = getwd(), limitsize = F)

#################################################
######        Supplementary Figures       #######
#################################################


# Missing data plots
for (m in c(0.1, 0.3, 0.5)) {
  b <- 1.5 * ((1/m) - 1)
  assign(paste0("miss.", m), rbeta(100000, shape1 = 1.5, shape2 = b))
}
miss <- data.frame(cbind(miss.0.1, miss.0.3, miss.0.5))
missing.perc.plot <- ggplot() + 
  geom_freqpoly(data = miss, mapping = aes(x=miss.0.1), bins = 25, color = "#c7e9b4", size = 1.5) +
  #geom_freqpoly(data = miss, mapping = aes(x=miss.0.2), bins = 25, color = "#7fcdbb", size = 1.5) +
  geom_freqpoly(data = miss, mapping = aes(x=miss.0.3), bins = 25, color = "#41b6c4", size = 1.5) +
  #geom_freqpoly(data = miss, mapping = aes(x=miss.0.4), bins = 25, color = "#2c7fb8", size = 1.5) +
  geom_freqpoly(data = miss, mapping = aes(x=miss.0.5), bins = 25, color = "#253494", size = 1.5) +
  xlab("Percent missing data")+
  ylab(paste("Proportion of individuals"))+
  labs(title = paste0("Distribution of missing data across individuals")) +
  xlim(c(0,1)) +
  scale_y_continuous(labels = function(x)x/100000, limits = c(0,30000)) +
  theme_classic() + 
  annotate(geom="text", x=0.32, y=22000, label="Overall missing data in dataset:", color="black") +
  annotate(geom="text", x=0.13, y=20000, label="0.1", color="black") +
  #annotate(geom="text", x=0.22, y=11500, label="0.2", color="black") +
  annotate(geom="text", x=0.35, y=8000, label="0.3", color="black") +
  #annotate(geom="text", x=0.45, y=7000, label="0.4", color="black") +
  annotate(geom="text", x=0.60, y=6200, label="0.5", color="black")
ggplot2::ggsave(filename = "missing.perc.plot.png", missing.perc.plot, width = 6, height = 4, path = getwd())




# Make plots for gifs

# triangle
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")
for (model in c("IBD", "SC", "TZ")) {
  if (model == "IBD") {
    full.model <- "IBD"
    new.model <- "IBD"
  }
  if (model == "SC") {
    full.model <- "SecondaryContact"
    new.model <- "ND"
  }
  if (model == "TZ") {
    full.model <- "TensionZone"
    new.model <- "GBC"
  }
  for (g in seq(from = 0, to = 19000, by =200)) {
    temp <- plotsForGifs(file = paste0(full.model,"/gen.",g,".vcf"), difference = 0.5, p1 = 0, p2 = 20, tri = T)
    ggplot2::ggsave(filename = paste0("Figures/AllTriangles/",new.model,".",g,".tri.png"), temp, 
                    width = 6, height = 4, path = getwd(), limitsize = F)
  }

}

# structure
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")
for (model in c("IBD", "SC", "TZ")) {
  if (model == "IBD") {
    full.model <- "IBD"
    new.model <- "IBD"
  }
  if (model == "SC") {
    full.model <- "SecondaryContact"
    new.model <- "ND"
  }
  if (model == "TZ") {
    full.model <- "TensionZone"
    new.model <- "GBC"
  }
  for (g in seq(from = 0, to = 19000, by = 200)) {
    dir <- paste0(full.model,"/str.output/01.output.",model,".gen.",g,".200K.K2/")
    files <- list.files(dir)
    files <- files[grepl("*f$", files)]  
    temp <- StructureImage(dir = dir, files = files, model.name = model, gen = g) 
    temp <- temp + labs(title = paste0("gen: ", g)) + theme(plot.title = element_text(size = 18, hjust = 0.5))
    ggplot2::ggsave(filename = paste0("Figures/AllStructures/",new.model,".",g,".str.png"), temp, 
                    width = 12, height = 4, path = getwd(), limitsize = F)
  }
}

# Make gifs

# triangle
IBD_triangles <- list.files("Figures/AllTriangles/", pattern = "IBD.\\d*.tri.png")
ND_triangles <- list.files("Figures/AllTriangles/", pattern = "ND.\\d*.tri.png")
GBC_triangles <- list.files("Figures/AllTriangles/", pattern = "GBC.\\d*.tri.png")

IBD_triangles <- IBD_triangles[order(nchar(IBD_triangles), IBD_triangles)]
ND_triangles <- ND_triangles[order(nchar(ND_triangles), ND_triangles)]
GBC_triangles <- GBC_triangles[order(nchar(GBC_triangles), GBC_triangles)]

IBD_triangles <- c(rep(IBD_triangles[1],3), IBD_triangles)
ND_triangles <- c(rep(ND_triangles[1],3), ND_triangles)
GBC_triangles <- c(rep(GBC_triangles[1],3), GBC_triangles)

gifski(paste0("Figures/AllTriangles/", IBD_triangles), gif_file ="Figures/IBD_triangles.gif", delay = 0.5, width = 600, height = 400)
gifski(paste0("Figures/AllTriangles/", ND_triangles), gif_file ="Figures/ND_triangles.gif", delay = 0.5, width = 600, height = 400)
gifski(paste0("Figures/AllTriangles/", GBC_triangles), gif_file ="Figures/GBC_triangles.gif", delay = 0.5, width = 600, height = 400)


# structure
IBD_structures <- list.files("Figures/AllStructures/", pattern = "IBD.\\d*.str.png")
ND_structures <- list.files("Figures/AllStructures/", pattern = "ND.\\d*.str.png")
GBC_structures <- list.files("Figures/AllStructures/", pattern = "GBC.\\d*.str.png")

IBD_structures <- IBD_structures[order(nchar(IBD_structures), IBD_structures)]
ND_structures <- ND_structures[order(nchar(ND_structures), ND_structures)]
GBC_structures <- GBC_structures[order(nchar(GBC_structures), GBC_structures)]

IBD_structures <- c(rep(IBD_structures[1],3), IBD_structures)
ND_structures <- c(rep(ND_structures[1],3), ND_structures)
GBC_structures <- c(rep(GBC_structures[1],3), GBC_structures)

gifski(paste0("Figures/AllStructures/", IBD_structures), gif_file ="Figures/IBD_structures.gif", delay = 0.5, width = 1200, height = 400)
gifski(paste0("Figures/AllStructures/", ND_structures), gif_file ="Figures/ND_structures.gif", delay = 0.5, width = 1200, height = 400)
gifski(paste0("Figures/AllStructures/", GBC_structures), gif_file ="Figures/GBC_structures.gif", delay = 0.5, width = 1200, height = 400)





# Cover image
gen6 <- read_table("triangle_plot_possible_space_math/gen6.txt", col_names = T)

cover.image <-  ggplot(gen6, aes(x=hi, y=het)) +
  geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
  geom_point(cex = 0.75, alpha = 1) +
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("Hybrid Index"))+
  ylab(paste("Heterozygosity"))+
  labs(title = "") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_classic() +
  geom_point(data = example.triangle, aes(x=hi, y=het, color = class), cex = 4) +
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  scale_color_manual(values = c("#e7d4e8", "#c2a5cf", "#9970ab", "#762a83",
                                "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837",  
                                "#bababa", "#bababa", "#40004b", "#00441b",
                                "#FFFFFF", "#FFFFFF", "#FFFFFF")) +
  theme_classic() +
  annotate("text", x=-0.075,y=-0.02,label="P1",color="#000000") +
  annotate("text", x=1.075,y=-0.02,label="P2",color="#000000") +
  annotate("text", x=0.5,y=1.075,label="F1",color="#000000") +
  annotate("text", x=0.5,y=0.425,label="F2 and later",color="#000000") +
  annotate("text", x=0.155,y=0.5,label="BC1",color="#000000") +
  annotate("text", x=0.845,y=0.5,label="BC1",color="#000000") +
  annotate("text", x=0.03,y=0.25,label="BC2",color="#000000") +
  annotate("text", x=-0.032,y=0.155,label="BC3",color="#000000") +
  annotate("text", x=-0.06075,y=0.0725,label="BC4",color="#000000") +
  annotate("text", x=0.97,y=0.25,label="BC2",color="#000000") +
  annotate("text", x=1.032,y=0.155,label="BC3",color="#000000") +
  annotate("text", x=1.06075,y=0.0725,label="BC4",color="#000000") +
  theme(axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16), axis.text=element_text(size=14)) +
  theme(legend.position = "none")

setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Cover_Image.pdf", cover.image, 
                  width = 5, height = 4, path = getwd(), limitsize = F)
  
ggplot2::ggsave(filename = "Cover_Image.png", cover.image, 
                width = 5, height = 4, path = getwd(), limitsize = F)









#    NOT RUN






# Missing data
for (f in c("IBD/gen.0.vcf", "IBD/gen.1000.vcf", "IBD/gen.10000.vcf")) {
  for (m in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
    generationPlotter(file = f, model = "IBD", missing.perc = m, 
                      missing.plot = T, write.plots = T, write.missing = T, write.str = T)
  }
}  

for (f in c("SecondaryContact/gen.0.vcf", "SecondaryContact/gen.1000.vcf", "SecondaryContact/gen.10000.vcf")) {
  for (m in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
    generationPlotter(file = f, model = "Secondary_Contact", missing.perc = m, 
                      missing.plot = T, write.plots = T, write.missing = T, write.str = T)
  }
} 

for (f in c("TensionZone/gen.0.vcf", "TensionZone/gen.1000.vcf", "TensionZone/gen.10000.vcf")) {
  for (m in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
    generationPlotter(file = f, model = "Tension_Zone", missing.perc = m, 
                      missing.plot = T, write.plots = T, write.missing = T, write.str = T)
  }
} 




########################################
####     Missing data Structure     ####
########################################


### Make IBD structure plots ###
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/IBD/str.miss.output")
# Make plot for every generation
dirs <- list.files()
for (dir in dirs) {
  setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/IBD/str.miss.output")
  generation <- sapply(strsplit(dir, "\\."),'[',5)
  missing <- paste0("0.", sapply(strsplit(dir, "\\."), '[', 8))
  files <- list.files(dir)
  files <- files[grepl("*f$", files)]
  saveStructureImage(files = files, model.name = "IBD", gen = generation, miss = missing, out_dir = "../str.miss.plots/")
}

### Make SC structure plots ###
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/SecondaryContact/str.miss.output")
# Make plot for every generation
dirs <- list.files()
for (dir in dirs) {
  setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/SecondaryContact/str.miss.output")
  generation <- sapply(strsplit(dir, "\\."),'[',5)
  missing <- paste0("0.", sapply(strsplit(dir, "\\."), '[', 8))
  files <- list.files(dir)
  files <- files[grepl("*f$", files)]
  saveStructureImage(files = files, model.name = "Secondary_Contact", gen = generation, miss = missing, out_dir = "../str.miss.plots/")
}

### Make TZ structure plots ###
setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/TensionZone/str.miss.output")
# Make plot for every generation
dirs <- list.files()
for (dir in dirs) {
  setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/TensionZone/str.miss.output")
  generation <- sapply(strsplit(dir, "\\."),'[',5)
  missing <- paste0("0.", sapply(strsplit(dir, "\\."), '[', 8))
  files <- list.files(dir)
  files <- files[grepl("*f$", files)]
  saveStructureImage(files = files, model.name = "Tension_Zone", gen = generation, miss = missing, out_dir = "../str.miss.plots/")
}


# Hybrid index triangle plots
generationPlotter(file = "SecondaryContact/gen.1000.vcf", model = "SC", fixed.diffs.triangle = T, p1 = 0, p2 = 20, difference = 0.5)




