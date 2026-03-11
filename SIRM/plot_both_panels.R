library(cowplot)

source("plot_empirical.R")
source("empirical/find_peak.R")

plots <- list(pie, panel2)
grobs <- lapply(plots, as_grob)
plot_widths <- lapply(grobs, function(x) {x$widths})
# Aligning the left margins of all plots
aligned_widths <- align_margin(plot_widths, "first")
# Aligning the right margins of all plots as well
aligned_widths <- align_margin(aligned_widths, "last")
# Setting the dimensions of plots to the aligned dimensions
for (i in seq_along(plots)) {
  grobs[[i]]$widths <- aligned_widths[[i]]
}
# Draw aligned plots
both_plots <- plot_grid(plotlist = grobs, ncol = 1, rel_heights = c(4, 1))

pdf("~/asr_writing/manuscript/figs/pie_prevalence_2week.pdf", width =8, height = 9)  
both_plots +  annotate("rect", xmin = 0, xmax = 1, ymin = .21, ymax = .2525,
                       fill = "white", alpha = 1)
dev.off()