library(ape) 
tree <- read.nexus("Ebola_glm_mascot.mcc.trees")

new_name <- unlist(lapply((strsplit(tree$tip.label, "|", fixed = TRUE)), function(lst) {lst[3]}))

tree$tip.label <- new_name

set.seed(1)
keep <- sample(new_name , size = 50, replace = FALSE)
subsample_tree <- keep.tip(tree, keep)
tree_label <- makeNodeLabel(subsample_tree, method = "number", prefix = "node")
tree_label$edge.length <- unlist(lapply(tree_label$edge.length,  function(x) {if (x < 0) {x <- 0.000001};  return(x)  }  ))

tree_label$edge.length <- tree_label$edge.length * 52 

write.csv(keep, "kept_samples.csv",  row.names = FALSE, quote = FALSE)
write.tree(tree_label, "out.1.tre")
