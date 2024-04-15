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

############################################################
###    Define function to make PCAs and Triangle plots   ###
############################################################

slopeTester <- function(file = NULL, model = NULL, 
                              downsample.pops = NA, inds.per.pop = 20,
                              missing.perc = NA, missing.plot = F,
                              fixed.diffs.triangle = F, p1 = NULL, p2 = NULL, difference = NULL,
                              plot = T) {
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
    
    if (fixed.diffs.triangle) {
      
      # get differences above allele frequency threshold
      vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = p1, p2 = p2, difference = difference)
      # get hybrid index
      hi <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = p1, p2 = p2)
      
      # number of differences above threshold
      d <- nrow(vcfR.diff@gt)
      
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
        geom_smooth(method = "lm", se = FALSE) +
        theme(legend.position = "none", axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))
      #theme(legend.position = c(0.13,0.75))
      
      print(diff.triangle)

    }
    
  }
  return(hi)
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





# Hybrid index triangle plots
for (sim in c("IBD", "SecondaryContact", "TensionZone")) {
  for (gen in c(0, 200, 1000, 10000, 19000)) {
    df <- slopeTester(file = paste0(sim, "/gen.", gen, ".vcf"), fixed.diffs.triangle = T, p1 = 0, p2 = 20, difference = 0.5)
    lm.test <- lm(het ~ hybrid.index, df)
    values <- summary(lm.test)$coefficients[2,]
    values <- c(sim, gen, values)
    if (!exists("sig.table")) {
      sig.table <- values
    } else {
      sig.table <- rbind(sig.table, values)
    }
    
  }
  sig.table <- data.frame(sig.table)
  colnames(sig.table) <- c("Sim", "Gen", "Slope", "SE", "t", "p.value")
}

sig.table$sim <- c(rep("IBD",5), rep("ND",5), rep("GBC", 5))
sig.table$Slope <- round(as.numeric(sig.table$Slope), 3)
sig.table$SE <- round(as.numeric(sig.table$SE), 3)
sig.table$t <- round(as.numeric(sig.table$t), 3)
sig.table$p.value <- signif(as.numeric(sig.table$p.value), 3)

write.table(sig.table, file = "Figures/triangle_plot_lm.txt", sep = "\t", col.names = T, row.names = F, append = F, quote = F)


