#install.packages("devtools")
#devtools::install_github("cmt2/RevGadgets")

library(RevGadgets)
library(ggplot2)

#####

file <- "~/asr_phyddle/liolaemus/output/1/lio_1_ase.tre"

tree <- readTrees(path=file)
plot <- plotTree(tree = tree, line_width = 0.5)
labs <- c("0" = "Andean", 
          "1" = "Lowland", 
          "2" = "Both")
geo_exam <- processAncStates(file, state_labels = labs)
pie <- plotAncStatesPie(t =geo_exam, cladogenetic = TRUE, 
                        tip_labels_states_offset = .05,
                        # Offset the tip labels to make room for tip pies
                        tip_labels_offset = .2, 
                        # Move tip pies right slightly 
                        tip_pie_nudge_x = .07,
                        state_transparency = 1.0,
                        # Change the size of node and tip pies  
                        tip_pie_size = 0.5,
                        node_pie_size = .5, 
                        ladderize = FALSE)
pdf("~/asr_writing/manuscript/figs/lio_geosse_Bayes.pdf")
print(pie)
dev.off()


file <- "~/asr_phyddle/liolaemus/phyddle/empirical/lio.1.est.tre"

tree <- readTrees(path=file)
plot <- plotTree(tree = tree, line_width = 0.5)
labs <- c("0" = "Andean", 
          "1" = "Lowland", 
          "2" = "Both")
#labs <- c("2" = "F", 
#          "3" = "C", 
#          "1" = "B")
geo_exam <- processAncStates(file, state_labels = labs)
pie <- plotAncStatesPie(t =geo_exam, cladogenetic = TRUE, 
                        tip_labels_states_offset = .05,
                        # Offset the tip labels to make room for tip pies
                        tip_labels_offset = .2, 
                        # Move tip pies right slightly 
                        tip_pie_nudge_x = .07,
                        # Change the size of node and tip pies  
                        tip_pie_size = 0.5,
                        node_pie_size = .5, 
                        state_transparency = 1.0,
                        ladderize = FALSE) 
pdf("~/asr_writing/manuscript/figs/lio_geosse_phyddle.pdf")
print(pie)
dev.off()

# #
# file <- "~/asr_phyddle/liolaemus/phyddle/empirical/lio.1.ann.nex"
# 
# tree <- readTrees(path=file)
# plot <- plotTree(tree = tree, line_width = 0.5)
# labs <- c("0" = "A->A,A",
#           "1" = "B->B,B",
#           "2" = "AB->A,B",
#           "3" = "AB->B,A",
#           "4" = "AB->AB,A",
#           "5" = "AB->A,AB",
#           "6" = "AB->AB,B",
#           "7" = "AB->B,AB")
# #labs <- c("2" = "F", 
# #          "3" = "C", 
# #          "1" = "B")
# geo_exam <- processAncStates(file, state_labels = labs)
# pie <- plotAncStatesPie(t =geo_exam, cladogenetic = FALSE, 
#                         tip_labels_states_offset = .05,
#                         # Offset the tip labels to make room for tip pies
#                         tip_labels_offset = .2, 
#                         # Move tip pies right slightly 
#                         tip_pie_nudge_x = .07,
#                         # Change the size of node and tip pies  
#                         tip_pie_size = 0.5,
#                         node_pie_size = .5, 
#                         ladderize = FALSE)
# pdf("~/asr_phyddle/liolaemus/phyddle.pdf")
# print(pie)
# dev.off()


ASR_map <- plotAncStatesMAP(t =geo_exam, 
                            cladogenetic = TRUE,
                            tip_labels_states_offset = .05,
                            # Offset the tip labels to make room for tip pies
                            # tip_labels_offset = .2, 
                            # Move tip pies right slightly 
                            #tip_pie_nudge_x = .07,
                            # Change the size of node and tip pies  
                            #tip_pie_size = 0.5,
                            #node_color= cols,
                            #node_pie_size = .5, 
                            ladderize = FALSE)
ASR_map
