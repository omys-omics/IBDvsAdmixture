setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")

library(readr)
library(patchwork)
library(ggpubr)
library(reshape)
library(ggplot2)

vignette(topic="model-comparison",package="conStruct")


xValPlotter <- function(model = NULL, gen = NULL) {
  if (m == "IBD") {
    newmodel <- "IBD"
  }
  if (m == "Secondary_Contact") {
    newmodel <- "ND"
  }
  if (m == "Tension_Zone") {
    newmodel <- "GBC"
  }
  
  # read in data
  sp.results <- as.matrix(
    read.table(paste0(model, ".gen.", gen, "_sp_xval_results.txt"),
               header = TRUE,
               stringsAsFactors = FALSE)
  )
  nsp.results <- as.matrix(
    read.table(paste0(model, ".gen.", gen, "_nsp_xval_results.txt"),
               header = TRUE,
               stringsAsFactors = FALSE)
  )
  
  # first, get the 95% confidence intervals for the spatial and nonspatial
  #   models over values of K (mean +/- 1.96 the standard error)
  
  sp.CIs <- data.frame(apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)}))
  colnames(sp.CIs) <- c("1", "2", "3")
  sp.CIs <- melt(sp.CIs, var="K")
  
  # combine data
  results <- cbind(sp.results, nsp.results)
  results <- data.frame(t(results))
  results$rep <- c(rep("spatial", 8), rep("non-spatial",8))
  colnames(results) <- c("1", "2", "3", "model")
  results <- melt(results, id="model", var="K")
  
  # plot results
  plot1 <- ggplot(results, aes(x=K, y=value, color = model)) +
    stat_summary(fun = "mean", geom = "point", cex = 2) +
    guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
    ggtitle(paste0(newmodel, " gen:", gen)) +
    xlab(paste("values of K"))+
    ylab(paste("predicitve accuracy"))+
    scale_color_manual(values = c("green", "blue")) +
    scale_y_continuous(limits = c(NA, 0)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    theme(legend.position = "none")
  
  plot2 <- ggplot() +
    stat_summary(data = results[results$model=="spatial",], aes(x=K, y=value, color = model), fun = "mean", geom = "point", cex = 2, color = "blue") +
    geom_line(data = sp.CIs, aes(x=K, y=value), color = "blue") +
    guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
    ggtitle("zoom in") +
    xlab(paste("values of K"))+
    ylab(paste("predicitve accuracy"))+
    scale_y_continuous(limits = c(NA, 0)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    theme(legend.position = "none")
  
  return(list(plot1, plot2))
}

setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/xVal")

for (m in c("IBD", "Secondary_Contact", "Tension_Zone")) {
  for (g in c(0, 200, 1000, 10000, 19000)) {
    if (m == "IBD") {
      n <- "IBD"
    }
    if (m == "Secondary_Contact") {
      n <- "ND"
    }
    if (m == "Tension_Zone") {
      n <- "GBC"
    }
    temp <- xValPlotter(model = m, gen = g)
    assign(paste0(n,".gen.",g,".plots"),temp)
  }
}

plot <- ggarrange(IBD.gen.0.plots[[1]] + IBD.gen.0.plots[[2]], ND.gen.0.plots[[1]] + ND.gen.0.plots[[2]], GBC.gen.0.plots[[1]] + GBC.gen.0.plots[[2]],
          IBD.gen.200.plots[[1]] + IBD.gen.200.plots[[2]], ND.gen.200.plots[[1]] + ND.gen.200.plots[[2]], GBC.gen.200.plots[[1]] + GBC.gen.200.plots[[2]],
          IBD.gen.1000.plots[[1]] + IBD.gen.1000.plots[[2]], ND.gen.1000.plots[[1]] + ND.gen.1000.plots[[2]], GBC.gen.1000.plots[[1]] + GBC.gen.1000.plots[[2]],
          IBD.gen.10000.plots[[1]] + IBD.gen.10000.plots[[2]], ND.gen.10000.plots[[1]] + ND.gen.10000.plots[[2]], GBC.gen.10000.plots[[1]] + GBC.gen.10000.plots[[2]],
          IBD.gen.19000.plots[[1]] + IBD.gen.19000.plots[[2]], ND.gen.19000.plots[[1]] + ND.gen.19000.plots[[2]], GBC.gen.19000.plots[[1]] + GBC.gen.19000.plots[[2]],
          
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
          font.label = list(size = 14),
          ncol = 3, nrow = 5)

setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "xVal.pdf", plot, 
                width = 11, height = 12, path = getwd(), limitsize = F)
ggplot2::ggsave(filename = "xVal.png", plot, 
                width = 11, height = 12, path = getwd(), limitsize = F)




