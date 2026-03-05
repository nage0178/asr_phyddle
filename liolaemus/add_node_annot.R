library(ape)

tree <- read.tree("data/lio.tre")

tree <- makeNodeLabel(tree, "n")
write.tree(tree, file = "phyddle/empirical/lio.1.tre")

