library(ggplot2)

#est_p <- read.csv(paste("../est.csv"), header = TRUE)
#est_p <- as.matrix(est_p)

#phyddle_mat <- matrix(NA, ncol = 8, nrow = 0) 
bayes_mat <- matrix(NA, ncol = 8, nrow = 0) 
truth_all <- c() 
for (i in 1:100) {
#for (i in 1:2500) {

	file2 <- paste("parseRb/node_index_", i, sep = "")
	if (! file.exists(file2)) {
		next
		print(paste("problem", i))
		q()
	}

	# Matches node names to indeces used in revbayes
	# The second column should have the right order
	mapping <- read.csv(paste("parseRb/node_index_", i, sep = ""))
	truth <- read.csv(paste("../simulate/sim.", i, ".anc_state.csv",  sep = ""), header = FALSE)

	# Estimate from rb
	# ANNA: Change to read in probs
	est_rb <- read.csv(paste("parseOutput/", i, "_PP.csv", sep = ""), header = TRUE)

       	order <- match(est_rb[,1], paste("ind_", mapping[,2], sep = ""))

	rb_reorder <- cbind(est_rb[order, ], mapping[,1]) 
	print(rb_reorder)
	truth <- truth[which(truth[,1] %in% rb_reorder[,10]), ]
	truth_all <- c(truth_all, truth[, 2])
	print(truth_all)

	bayes_mat <- rbind(bayes_mat, rb_reorder)
	#q()


}
write.csv(bayes_mat, "bayes_mat.csv")
bayes_mat <- read.csv("bayes_mat.csv")
rownames(bayes_mat) <- bayes_mat[,1]
bayes_mat <- bayes_mat[, -1]
head(bayes_mat)
print(dim(bayes_mat))

prob_true <- function(row){
  true_state <- row[ 9]
  row[true_state] < .05
}

#conf_wrong <- matrix(NA, nrow = dim(truth)[1], ncol= n_node)
tmp <- cbind(bayes_mat, truth_all)
tmp <- tmp[, -1]
tmp <- tmp[, -9]
head(tmp)
conf_wrong <- apply(tmp, 1, prob_true)

# double check you are comparing the right states
for (i in 1:8) {
  state <- which(truth_all == i)
  print(mean(conf_wrong[state]))
}
mean (conf_wrong)
