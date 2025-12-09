IDX <- 1

#write.table(mcmc[, -c(1,2)], "test2.txt", row.names = FALSE)
#q()

set.seed(1)

# Determines the state based on the parent/daughter triplet
# The numbering matches that used in phyddle
find_triple <- function(mat, parent, d_l, d_r) {
	if (mat[parent] == 0 && mat[d_l] == 0 && mat[d_r] == 0) {
		state <- 1
	} else if (mat[parent] == 1 && mat[d_l] == 1 && mat[d_r] == 1) {
		state <- 2
	} else if (mat[parent] == 2 && mat[d_l] == 0 && mat[d_r] == 1) {
		state <- 3
	} else if (mat[parent] == 2 && mat[d_l] == 1 && mat[d_r] == 0) {
		state <- 4
	} else if (mat[parent] == 2 && mat[d_l] == 2 && mat[d_r] == 0) {
		state <- 5
	} else if (mat[parent] == 2 && mat[d_l] == 0 && mat[d_r] == 2) {
		state <- 6
	} else if (mat[parent] == 2 && mat[d_l] == 2 && mat[d_r] == 1) {
		state <- 7
	} else if (mat[parent] == 2 && mat[d_l] == 1 && mat[d_r] == 2) {
		state <- 8
	}
	state
}

for (IDX in 1:2500) {

	file <- paste("output/1/", IDX, "_ase.tre", sep = "")

	if ( file.exists(file)) {
		print(IDX)
		mcmc <- read.table(paste("output/1/", IDX, "_states.log", sep = ""), header = TRUE)

		# Matrix of parent daughter relationships with revbayes indeces
		par_mat <- read.table(paste("output/1/", IDX, "_matrix.txt", sep = ""), header = TRUE)
		new_anc_states <- matrix(NA, nrow = dim(mcmc)[1], ncol = dim(par_mat)[1])
		
		# Names are based on the parent node index
		colnames(new_anc_states) <- par_mat[,1] 
		
		# For every internal node in the tree
		for (i in 1:dim(par_mat)[1]) {

			# Finds columns corresponding to start and end of parent branch
			par <- par_mat[i, 1]
			par_start <- 2 + ((par - 1) * 2) #+ 1 #This + 1 is if there are two replicates
			par_end <- par_start + 1
		
		
			# Finds columns corresponding to start the two daughter branchs
			d <- par_mat[i, 2:3]
			dau_1_start <- 2 + ((d[[1]] - 1) * 2) #+ 1 #If two replicates in same file need this
			dau_2_start <- 2 + ((d[[2]] - 1) * 2) #+ 1
		
			
			# Determines the state for each row in MCMC
			output <- apply(mcmc, 1, find_triple, parent = par_end, d_l =dau_1_start , d_r = dau_2_start)
			new_anc_states[, i] <- output
			colnames(new_anc_states)[i] <- paste("ind_", par, sep = "")
		}

		# Table with state for each row in MCMC for each internal node 
		write.table(new_anc_states, paste("parseOutput/", IDX, "_anc_state.txt", sep = ""), sep = "\t", row.names = FALSE)
		
		# Finds the point estimate based on which state is most frequent
		point_est <- matrix(NA, nrow = dim(new_anc_states)[2], ncol = 2)
		for ( i in 1:dim(new_anc_states)[2]) {
			sum_inf <- table(new_anc_states[,i])
			colNum <- which(sum_inf == max(sum_inf))

			point_est[i, 1] <- colnames(new_anc_states)[i]
			most_common <- as.numeric(names(sum_inf)[colNum])
			if (length(most_common) > 1) {
				point_est[i, 2] <- sample(most_common, 1)
			} else { 
				point_est[i, 2] <- most_common
			}
		}
		write.table(point_est, paste("parseOutput/", IDX, "point_est.txt", sep = ""), sep = ",", row.names = FALSE, col.names=FALSE, quote = FALSE)
	} else {
		print(paste("problem with ", IDX))
		q()
	}
}
