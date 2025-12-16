library(ggplot2)
library(ggpubr)
library(patchwork)

source("~/asr_phyddle/mk_binary/checkPerformance.R")
source("~/asr_phyddle/mk_binary/plot_node_accuracy.R")
source("~/asr_phyddle/mk_binary/fixedSizeInVarNet.R")

scaling <- 1
pdf("~/asr_writing/manuscript/figs/fix_var.pdf", height = 5 * scaling, width= 8 * scaling )
ggarrange(
          ggarrange(panelA, panelB, ncol = 2, labels = c("a", "b")), # Second row with box and dot plots
          panelC,                                                 # First row with scatter plot
          nrow = 2, 
          labels = c("", "c")                                        # Labels of the scatter plot
) 

dev.off()
