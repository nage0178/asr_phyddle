library(RevGadgets)
library(ggplot2)
library(ggpubr)

file <- paste("~/asr_phyddle/SIRM/empirical/out.1.est.tre", sep = "")
#file <- paste("~/asr_phyddle/SIRM_original_script/vhigh_sample/smaller_sample/migration/50_tip/empirical/out.1.est.tre", sep = "")


tree <- readTrees(path=file)
  #plot <- plotTree(tree = tree, line_width = 0.5)
labs <- c("0" = "0",
            "1" = "1",
            "2" = "2", 
            "3" = "3",
            "4" = "4")

cols <- c("darkred", "darkviolet", "darkgoldenrod1", "darkgreen", "blue1")
  
geo_exam <- processAncStates(file , state_labels = labs)
pie <- plotAncStatesPie(t =geo_exam, 
                          tip_labels_states_offset = .05,
                          # Offset the tip labels to make room for tip pies
                          # tip_labels_offset = .2, 
                          # Move tip pies right slightly 
                          tip_pie_nudge_x = .07,
                          # Change the size of node and tip pies  
                          tip_pie_size = 0.5,
                          pie_colors= cols,
                          node_pie_size = .75, 
                        state_transparency = 1.0, 
                          ladderize = FALSE)

pdf("~/asr_writing/manuscript/figs/Ebola_pie.pdf", width =8, height = 6)  
print(pie)
dev.off()


pdf("~/asr_writing/manuscript/figs/Ebola_pie_alone.pdf", width =8, height = 6)  
pie + theme(legend.position = "none")
dev.off()
system("pdfcrop ~/asr_writing/manuscript/figs/Ebola_pie_alone.pdf")


legend_pie <- cowplot::get_legend(pie)
leg_pie <- as_ggplot(legend_pie)

pdf("~/asr_writing/manuscript/figs/Ebola_leg_pie_alone.pdf", width =2, height = 6.75)  
leg_pie
dev.off()
  
ASR_map <- plotAncStatesMAP(t =geo_exam, 
                          tip_labels_states_offset = .05,
                          # Offset the tip labels to make room for tip pies
                          # tip_labels_offset = .2, 
                          # Move tip pies right slightly 
                          #tip_pie_nudge_x = .07,
                          # Change the size of node and tip pies  
                          #tip_pie_size = 0.5,
                          node_color= cols,
                          #node_pie_size = .5, 
                          state_transparency = 1.0, 
                          ladderize = FALSE)



library(cowplot)
map_alone <- ASR_map + theme(legend.position = "none")
legend1 <- cowplot::get_legend(ASR_map)
leg1 <- as_ggplot(legend1)

pdf("~/asr_writing/manuscript/figs/Ebola_leg_alone.pdf", width =2, height = 6.75)  
leg1
dev.off()

pdf("~/asr_writing/manuscript/figs/Ebola_map_alone.pdf", width =8, height = 6.75)  
map_alone 
dev.off()
system("pdfcrop ~/asr_writing/manuscript/figs/Ebola_map_alone.pdf")
