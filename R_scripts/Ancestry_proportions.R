setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD")

library(readr)
library(patchwork)
library(ggpubr)
library(reshape)

# Read in and format ancestry proportion data
IBD_anc_pro <- read_delim("Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/IBD.anc.pro_ORIGINAL.txt", delim = ":", escape_double = FALSE, trim_ws = TRUE, col_names = F)
IBD_anc_pro <- IBD_anc_pro[seq(2,84, by = 2),]
IBD_anc_pro <- cbind(IBD_anc_pro[1:21,2],IBD_anc_pro[22:42,2])
colnames(IBD_anc_pro) <- c("IBD_p0_anc", "IBD_p20_anc")
rownames(IBD_anc_pro) <- c(0:20)

ND_anc_pro <- read_delim("Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/secondary_contact.anc.pro_ORIGINAL.txt", delim = ":", escape_double = FALSE, trim_ws = TRUE, col_names = F)
ND_anc_pro <- ND_anc_pro[seq(2,84, by = 2),]
ND_anc_pro <- cbind(ND_anc_pro[1:21,2],ND_anc_pro[22:42,2])
colnames(ND_anc_pro) <- c("ND_p0_anc", "ND_p20_anc")
rownames(ND_anc_pro) <- c(0:20)

GBC_anc_pro <- read_delim("Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/tension_zone.anc.pro_ORIGINAL.txt", delim = ":", escape_double = FALSE, trim_ws = TRUE, col_names = F)
GBC_anc_pro <- GBC_anc_pro[seq(2,84, by = 2),]
GBC_anc_pro <- cbind(GBC_anc_pro[1:21,2],GBC_anc_pro[22:42,2])
colnames(GBC_anc_pro) <- c("GBC_p0_anc", "GBC_p20_anc")
rownames(GBC_anc_pro) <- c(0:20)

anc_pro <- cbind(0:20, IBD_anc_pro, ND_anc_pro, GBC_anc_pro)
colnames(anc_pro) <- c("pop", colnames(anc_pro)[-1])

anc_pro <- melt(anc_pro, id.vars="pop" ,var="anc_type")
anc_pro$value <- as.numeric(anc_pro$value)
anc_pro$sim <- c(rep("Isolation-By-Distance", 42), rep("Neutral Diffusion", 42), rep("Geographically-Bounded Contact", 42))
anc_pro$anc_type <- rep(c(rep("p0", 21), rep("p20", 21)), 3)

plot <- ggplot(anc_pro, aes(x=pop, y=value, facet = sim, color = anc_type)) +
  facet_grid(. ~ factor(sim, levels = c("Isolation-By-Distance", "Neutral Diffusion", "Geographically-Bounded Contact"))) +
  geom_line(cex = 1, alpha=1)+
  geom_point(cex = 1, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("Population"))+
  ylab(paste("% Ancestry from Parental Pops"))+
  scale_color_manual("Ancestry Type", values = c("#313695", "#a50026")) +
  theme_classic() +
  theme(legend.position = "right", axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12))


setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/SecondaryContact_vs_IBD/Figures")
ggplot2::ggsave(filename = "Ancestry_Proportions_gen19000.pdf", plot, 
                width = 8, height = 3, path = getwd(), limitsize = F)
ggplot2::ggsave(filename = "Ancestry_Proportions_gen19000.png", plot, 
                width = 8, height = 3, path = getwd(), limitsize = F)


