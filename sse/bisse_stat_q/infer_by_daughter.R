library(phytools)
library(ggplot2)

# This is for a naive estimator that estimates the state at a node
# as the most common daughter state"
inference <- matrix(NA, nrow = 2500 * 49, ncol = 4)
row <- 1
set.seed(1)
dir <- "simulate/"

# For the test dataset
for (j in 1:2500) {
  tree <- read.newick(paste(dir, "out.", j,  ".tre", sep = ""))
  dat <- read.csv(paste(dir, "out." , j, ".dat.csv", sep = ""))
  ages <- read.csv(paste(dir, "out.", j, ".age_prop.csv", sep = ""))
  truth <- read.csv(paste(dir, "out.", j, ".anc_state.csv", sep = ""), header = FALSE)
    
  
 # For all the internal nodes (they start at n + 1)
 for (i in 51:99) {
    all_descend <- getDescendants(tree, i)
    tip_descend <- all_descend[which(all_descend < 51)]
    tip_nodes <- as.numeric(tree$tip.label[tip_descend])

    # Infer the node state based on the tip states alone
    if ( mean(dat[tip_nodes, 2]) > .5 ) {
      inf <- 1
    } else if (mean(dat[tip_nodes, 2]) < .5 ) {
      inf <- 0
    } else {
      inf <- sample(0:1, 1)
    }
    
    # Find the node name
    node <- tree$node.label[i-50]
    # Find the node height based on the name
    height <- ages[node]
    # Check if inference is correct
    correct <- (inf == truth[which(truth[, 1] == node), 2])

    inference[row, 1] <- j
    inference[row, 2] <- i
    inference[row, 3] <- correct
    inference[row, 4] <- as.numeric(height)
    
    row <- row + 1
  }
  
}

mean_correct <- matrix(0, nrow = 11, ncol = 2)
colnames(mean_correct) <- c("avg_age", "avg_correct")

for (i in 0:10) {
  index <- which(i == floor(inference[, 4] * 10))
  
  mean_correct[i+1, 1] <- mean(inference[index,4])
  mean_correct[i+1, 2] <- mean(inference[index,3])

}

df <- as.data.frame(mean_correct)
write.csv(df, "prop_tips_clade_inference_correct.csv")
#ggplot(df, aes(x = avg_age, y = avg_correct)) + geom_point() + geom_line()

