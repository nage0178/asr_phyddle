for (i in 1:2500) {
	mcmc <- read.table(paste("parseOutput/", i, "_anc_state.txt", sep = ""), sep = "\t", header = TRUE)
	df <- as.data.frame(mcmc)
	for (j in 1:ncol(mcmc)) {
		df[, j] <- factor(mcmc[, j], levels = c("1","2","3","4","5","6","7","8"))
	}


	freq <- sapply(df, function(x) xtabs(~x,  drop.unused.levels = FALSE), simplify = FALSE)
	freq <- matrix(unlist(freq), byrow = TRUE, nrow = ncol(mcmc)) / nrow(mcmc)
	colnames(freq) <- c(1:8)
	rownames(freq) <- colnames(mcmc)
	#print(freq)
	write.csv(freq, paste("parseOutput/", i, "_PP.csv", sep = ""))
}
