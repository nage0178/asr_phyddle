library(RevGadgets)
library(ggplot2)
library(ggtree)
library(ggpubr)

for (i in 5:1) {
#file <- paste("~/asr_phyddle/SIRM_original_script/vhigh_sample/smaller_sample/migration/50_tip/empirical/train_", i, "_out.3.est.tre", sep = "")

file <- paste("~/asr_phyddle/SIRM/empirical/train_", i, "_out.1.est.tre", sep = "")

labs <- c("0" = "0",
            "1" = "1",
            "2" = "2", 
            "3" = "3",
            "4" = "4")

cols <- c("darkred", "darkviolet", "darkgoldenrod1", "darkgreen", "blue1")

  
geo_exam <- processAncStates(file , state_labels = labs)

geo_exam@phylo$edge.length <- geo_exam@phylo$edge.length /52
pie <- plotAncStatesPie(t =geo_exam, 
                          #tip_labels_states_offset = .02,
                          # Offset the tip labels to make room for tip pies
                          tip_labels_offset = .002, 
                          tip_labels_size = 3,
                          tip_labels = FALSE,
                          # Move tip pies right slightly 
                          tip_pie_nudge_x = .00,
                          # Change the size of node and tip pies  
                          tip_pie_size = 0.5,
                          pie_colors= cols,
                          node_pie_size = .75, 
                          state_transparency = 1.0, 
                          ladderize = FALSE, 
                          mrsd="2014-12-08") +theme_tree2() + ggplot2::scale_x_continuous(c(-10,10)) 
if (i == 1) {
  pie <- pie + xlim_tree(2015 + 1/12)
}


pdf(paste("~/asr_writing/manuscript/figs/prev_train_", i, "_Ebola_pie.pdf", sep = ""), width =8, height = 6)  
print(pie)
dev.off()
cmd <- paste("pdfcrop --margins '0 0 0 -10' ~/asr_writing/manuscript/figs/prev_train_", i, "_Ebola_pie.pdf", sep = "")
system(cmd)

}

