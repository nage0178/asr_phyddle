library(ggplot2)

est_p <- read.csv(paste("../est.csv"), header = TRUE)
est_p <- as.matrix(est_p)

phyddle_mat <- matrix(NA, ncol = 8, nrow = 0) 
bayes_mat <- matrix(NA, ncol = 8, nrow = 0) 
for (i in 1:2500) {

	file <- paste("parseOutput/", i, "point_est.txt", sep = "")
	file2 <- paste("parseRb/node_index_", i, sep = "")
	if (! file.exists(file) || ! file.exists(file2)) {
		next
		print(paste("problem", i))
		q()
	}

	# Matches node names to indeces used in revbayes
	# The second column should have the right order
	mapping <- read.csv(paste("parseRb/node_index_", i, sep = ""))

	# Estimate from rb
	# ANNA: Change to read in probs
	est_rb <- read.csv(paste("parseOutput/", i, "_PP.csv", sep = ""), header = TRUE)

       	order <- match(est_rb[,1], paste("ind_", mapping[,2], sep = ""))

	rb_reorder <- est_rb[order, -1 ]
	bayes_mat <- rbind(bayes_mat, rb_reorder)

	row <- est_p[i, 2:ncol(est_p)]
	est_p_one <- matrix(row, byrow = TRUE, nrow = 49, ncol = 8)
	phyddle_mat <- rbind(phyddle_mat, est_p_one)

	#break

}
#dim(bayes_mat)
#dim(phyddle_mat)
method_diff <- phyddle_mat - bayes_mat 
write.csv(method_diff, "difference_phy-bayes.csv")
#colnames(method_diff)

library(tidyr)

long <- method_diff %>% 
  pivot_longer(
    cols = `X1`:`X8`, 
    names_to = "state",
    values_to = "diff"
  )
long

head(long)
df <- as.data.frame(long)
ggplot(df, aes(x = diff)) + geom_histogram() + facet_wrap(~state)

