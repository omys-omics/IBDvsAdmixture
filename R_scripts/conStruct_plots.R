setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/conStructPlots")

library(readr)
library(patchwork)
library(ggpubr)
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
    ggplot2::ylab("Anc. Pro") +
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

plotter <- function(model = NULL, gen = NULL) {
  if (model == "IBD") {
    n <- "IBD"
  }
  if (model == "Secondary_Contact") {
    n <- "ND"
  }
  if (model == "Tension_Zone") {
    n <- "GBC"
  }
  
  load(paste0(model, ".gen.", gen, "_conStruct.results.Robj"))
  
  admix.props <- conStruct.results$chain_1$MAP$admix.proportions

  cols <- c("#313695", "#a50026")
  
  qmat <- data.frame(cbind(0:20, rep(0,21), 0:20, admix.props))
  colnames(qmat) <- c("id", "pct.miss", "orig.pop", "Group.1", "Group.2")
  # make pop 0 blue
  if (qmat[1,"Group.1"] < 0.5) {
    qmat$Group.1 <- 1 - qmat$Group.1
    qmat$Group.2 <- 1 - qmat$Group.2
  }
  conStructPlot <- betterStructurePlot(qmat, label.pops = F, horiz = F, legend.position = "none", col = cols, sort.probs = F, type = "bar")
  
  #conStructPlot <- conStructPlot + labs(title = paste0(model, ", gen: ", gen))
  
  return(conStructPlot)
  
}

#######################################
####      End define functions     ####
#######################################


# how to visualize the output of a conStruct model
vignette(topic="visualize-results",package="conStruct")


setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/conStructPlots")

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
    temp <- plotter(model = m, gen = g)
    assign(paste0(n,".gen.",g,".plot"),temp)
  }
}



plot <- ggarrange(IBD.gen.0.plot, ND.gen.0.plot, GBC.gen.0.plot,
          IBD.gen.200.plot, ND.gen.200.plot, GBC.gen.200.plot,
          IBD.gen.1000.plot , ND.gen.1000.plot, GBC.gen.1000.plot,
          IBD.gen.10000.plot , ND.gen.10000.plot, GBC.gen.10000.plot,
          IBD.gen.19000.plot , ND.gen.19000.plot, GBC.gen.19000.plot,
          
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
          font.label = list(size = 14),
          ncol = 3, nrow = 5,)

newplot <- annotate_figure(plot, 
                top = "Isolation-By-Distance                         Neutral Diffusion                  Geographically-Bounded Contact")


setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "all_conStruct.pdf", newplot, 
                width = 10, height = 6, path = getwd(), limitsize = F)




