library(viridis)
library(ggplot2)

sym1 <- read.csv("~/asr_phyddle/4_tip/3_bin/sym_train1.csv")
sym2 <- read.csv("~/asr_phyddle/4_tip/3_bin/sym_train2.csv")
#reorder <- c(1,8, 2,7, 5, 4, 3,6)
reorder <- c(1,3, 7, 6, 5, 8, 4 ,2)
sym1$order <- reorder
sym2$order <- reorder

library(data.table)
long_sym1 <- melt(setDT(sym1), id.vars = c("order"), variable.name = "method")
long_sym2 <- melt(setDT(sym2), id.vars = c("order"), variable.name = "method")

long_sym1$train <- "1"
long_sym2$train <- "2"
sym <- rbind(long_sym1, long_sym2)

sym$order <- as.factor(sym$order)
sym$method_train <- interaction(sym$train, sym$method)
#pdf("~/asr_writing/manuscript/figs/.pdf", width =4, height = 3)
# plotting the matrix
plot2 <-ggplot(sym, aes(y =method_train, x =  order, fill = value)) +
  geom_tile() + geom_text(aes(label = round(value, digits = 2)), col = "red") + scale_fill_viridis() + labs(x = "tree", y = "method, training", fill = "accuracy") +
  scale_y_discrete(labels = c("2.joint" = "joint, 2", 
                              "1.joint" = "joint, 1", 
                              "2.single.node" = "single node, 2", 
                              "1.single.node" = "single node, 1" , 
                              "2.marginal" = "marginal, 2" , 
                              "1.marginal" = "marginal, 1")) + 
  #geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'red', size = 8) +
  theme_light()
#dev.off()


  
library(RevGadgets)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree1.nex")
tree1 <- plotTree(tree, node_labels = "state", 
         node_labels_offset = 0.050,   line_width = 0.5,
         # italicize tip labels 
         tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree2.nex")
tree2 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree3.nex")
tree3 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree4.nex")
tree4 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree5.nex")
tree5 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree6.nex")
tree6 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree7.nex")
tree7 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/tmp/tree8.nex")
tree8 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
plot1 <- ggarrange(tree1, tree2, tree5, tree3, tree8, tree7, tree4, tree6, nrow = 2, ncol = 4, 
          labels = c("1", "3", "5", "7", "2", "4", "6", "8"), font.label = list(size = 8), label.x = .4)
pdf("~/asr_writing/manuscript/figs/4_tip_heatmap.pdf", width =10, height = 3)
ggarrange(plot1, plot2, ncol = 2, nrow = 1, labels = c("a", "b"), widths = c(.4, .6))
dev.off()

png("~/asr_writing/manuscript/figs/4_tip_heatmap.png", width =10, height = 3, units = "in", res = 500)
ggarrange(plot1, plot2, ncol = 2, nrow = 1, labels = c("", ""), widths = c(.4, .6))
dev.off()

tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree1.nex")
asy_1 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree2.nex")
asy_2 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree3.nex")
asy_3 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree4.nex")
asy_4 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree5.nex")
asy_5 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree6.nex")
asy_6 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree7.nex")
asy_7 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)
tree <- readTrees("~/asr_phyddle/4_tip/3_bin/trees/tree8.nex")
asy_8 <- plotTree(tree, node_labels = "state", 
                  node_labels_offset = 0.050,   line_width = 0.5,
                  # italicize tip labels 
                  tip_labels_italics = TRUE)

plot8 <- ggarrange(asy_1, asy_5, asy_2, asy_3, asy_8, asy_4, asy_7, asy_6, nrow = 2, ncol = 4, 
                   labels = c("1", "3", "5", "7", "2", "4", "6", "8"), font.label = list(size = 8), label.x = .4)
asym1 <- read.csv("~/asr_phyddle/4_tip/3_bin/asym_train1.csv")
asym2 <- read.csv("~/asr_phyddle/4_tip/3_bin/asym_train2.csv")
#reorder <- c(1,8, 2,7, 5, 4, 3,6)
reorder <- c(1, 5, 7, 4, 3, 8, 6 ,2)
asym1$order <- reorder
asym2$order <- reorder

library(data.table)
long_asym1 <- melt(setDT(asym1), id.vars = c("order"), variable.name = "method")
long_asym2 <- melt(setDT(asym2), id.vars = c("order"), variable.name = "method")

long_asym1$train <- "1"
long_asym2$train <- "2"
asym <- rbind(long_asym1, long_asym2)

asym$order <- as.factor(asym$order)
asym$method_train <- interaction(asym$train, asym$method)
# plotting the matrix
plot2 <- ggplot(asym, aes(y =method_train, x =  order, fill = value)) +
  geom_tile() + geom_text(aes(label = round(value, digits = 2)), col = "red") + scale_fill_viridis() + labs(x = "tree", y = "method, training", fill = "accuracy") +
  scale_y_discrete(labels = c("2.joint" = "joint, 2", 
                              "1.joint" = "joint, 1", 
                              "2.single.node" = "single node, 2", 
                              "1.single.node" = "single node, 1" , 
                              "2.marginal" = "marginal, 2" , 
                              "1.marginal" = "marginal, 1")) + 
  #geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'red', size = 8) +
  theme_light()
plot2


pdf("~/asr_writing/manuscript/figs/4_tip_heatmap_asym.pdf", width =10, height = 3)
ggarrange(plot8, plot2, ncol = 2, nrow = 1, labels = c("a", "b"), widths = c(.4, .6))
dev.off()

png("~/asr_writing/manuscript/figs/4_tip_heatmap_asym.png", width =10, height = 3, units = "in", res = 500)
ggarrange(plot8, plot2, ncol = 2, nrow = 1, labels = c("a", "b"), widths = c(.4, .6))
dev.off()

dev.off()

