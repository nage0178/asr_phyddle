comb_All <- matrix(NA, ncol = 4, nrow = 0)
#for (i in 1:100) {
for (i in 1:2500) {

	file <- paste("parseOutput/", i, "point_est.txt", sep = "")
	file2 <- paste("parseRb/node_index_", i, sep = "")
	if (! file.exists(file) || ! file.exists(file2)) {
		next
		print(paste("problem", i))
		q()
	} 

	# Matches node names to indeces used in revbayes
	mapping <- read.csv(paste("parseRb/node_index_", i, sep = ""))

	# Estimate from rb
	est <- read.csv(paste("parseOutput/", i, "point_est.txt", sep = ""), header = FALSE)

	# True state
	truth <- read.csv(paste("../simulate/sim.", i, ".anc_state.csv", sep = ""), header = FALSE)

	# Finds how to match the two orders
       	order <- match(est[,1], paste("ind_", mapping[,2], sep = ""))

	reorder <- cbind(mapping[order, ], est[,2])

	# Reorders truth 
	truth_order <- match(reorder[,1], truth[,1])
	all <- cbind(reorder, truth[truth_order,])

	comb_All <- rbind(comb_All, all)
	
}
colnames(comb_All)[3] <- "est"
colnames(comb_All)[5] <- "true"
print(comb_All)

df <- comb_All[, c(3,5)]
df$true <- as.factor(df$true)
df$est <- as.factor(df$est)
table(df$est)
table(df$true)
library(caret)
print("proportion correct")
print(mean(df$est == df$true))
cm <- confusionMatrix(df$est, df$true)
cm
cm_d <- as.data.frame(cm$table)

write.csv(cm_d, "bayes_confusion.csv", row.names = FALSE)
