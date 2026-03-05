library(ape)
tree <- read.tree("empirical/lio.1.form.tre")
#print(tree)
tips <- tree$Nnode +1
parent_daughter <- matrix(NA, nrow = tips - 1, ncol = 3)
for (i in (tips +1):(2*tips - 1)) {
  rows <- which(tree$edge[, 1] == i)
  names <- c(NA, NA)

  for (j in 1:2) {
    # Daughter is internal node
    index <- tree$edge[rows[j], 2]
    if (index > tips) {
      names[j] <- tree$node.label[index - tips]
      # Daugther is tip
    } else {
      names[j] <- tree$tip.label[index]
    }
  }
  parent_daughter[i - tips, ] <- c(tree$node.label[i - tips], names)
}

colnames(parent_daughter) <- c("parent", "left", "right")
write.csv(parent_daughter, "empirical/parent_daughter.csv", quote = FALSE, row.names = FALSE )
