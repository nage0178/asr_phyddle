library(treebalance)
library(ape)
library(rhdf5)
reps <- 25000 * 2
asymm <- rep(NA, reps)

for (i in 1:reps) {
	j <- i + 150000
	tree <- read.tree(paste("accuracy/out.", format(j, scientific = FALSE), ".tre", sep = ""))
	asymm[i] <- collessI(tree) == 3
}

est <- read.csv("estimate/out.empirical_est.labels_cat.csv")
H5_truth <- H5Fopen("format_accuracy/format/out.train.hdf5")
truth <- t(H5_truth$labels)

asr_0 <- truth[, 1] == 1
asr_1 <- truth[, 2] == 1
asr_2 <- truth[, 3] == 1

states <- matrix(NA, nrow = reps, ncol = 8 )
# Determine all three states
states[,1] <- (!asr_0) * (!asr_1) * (!asr_2)  
states[,2] <- (!asr_0) * (!asr_1) * ( asr_2)  
states[,3] <- (!asr_0) * ( asr_1) * (!asr_2)  
states[,4] <- (!asr_0) * ( asr_1) * ( asr_2)  
states[,5] <- ( asr_0) * (!asr_1) * (!asr_2)  
states[,6] <- ( asr_0) * (!asr_1) * ( asr_2)  
states[,7] <- ( asr_0) * ( asr_1) * (!asr_2)  
states[,8] <- ( asr_0) * ( asr_1) * ( asr_2)  

asymmMat <- matrix(NA, nrow = 8, ncol = 4)
symmMat <- matrix(NA, nrow = 8, ncol = 4)
print("done reading trees")
ml_est <- rep(NA, reps/2)
for (i in 1:(reps/2)) {
	j <- i + 150000
	true_state <- read.csv(paste("accuracy/out.", format(j, scientific = FALSE), ".anc_state.csv", sep = ""), header= FALSE)
	infer_state <- read.csv(paste("accuracy/out.", format(j, scientific = FALSE), ".asr.csv", sep = ""), header = FALSE)
	ml_est[i] <- all(true_state[, 2] == infer_state[,2])
}

print("ML average")
mean(ml_est)
for (j in 0:1) {
	# Asymmetric is second
	rows1 <- which((asymm*1) == j) 
	for (i in 1:8) {
		rows2 <- which (states[, i] == 1)
		rows <- intersect(rows1, rows2)
		rows <- which(rows < 175001)
		if (j == 0 ) {
			symmMat[i, 4] <- mean(ml_est[rows])
		} else {
			asymmMat[i, 4] <- mean(ml_est[rows])
		}
	}
	print("")
}

point_est <- matrix(NA, nrow = reps, ncol = 3)

nstates <- 3
for (i in c(1:nstates)) {
  colZero <- 2+(i-1)*2
  colOne <- 3+(i-1)*2
  
  # There shouldn't be any NAs since the trees are the same size for a dataset 
  point_est[, i] <- (est[,colOne] == apply(est[,colZero:colOne], 1, max)) * 1

}
correct <- apply(point_est == truth[,1:3], 1, all)
length(correct)
print("correct 3")
mean(correct)
print("by states")

for (j in 0:1) {
	# Asymmetric is second
	rows1 <- which((asymm*1) == j) 
	for (i in 1:8) {
		rows2<- which (states[, i] == 1)
		rows <- intersect(rows1, rows2)
		if (j == 0 ) {
			symmMat[i, 1] <- mean(correct[rows])
		} else {
			asymmMat[i, 1] <- mean(correct[rows])
		}
	}
}
# by node

H5Fclose(H5_truth)
print("starting one")
est <- read.csv("../one/estimate/est2.empirical_est.labels_cat.csv")
point_estOne <- matrix(NA, nrow = dim(truth)[1], ncol = 3)
point_estOne_vec <- (est[,3] ==  apply(est[,2:3], 1, max)) * 1

for (i in 1:nstates) {
  rows <- seq(i, reps * 3, by= 3)
  point_estOne[, i] <- point_estOne_vec[rows]
}

print("by node match")
mean(point_estOne == truth[,1:3])
correct <- apply(point_estOne == truth[,1:3], 1, all)
print("correct one")
mean(correct)
print("by states")

for (j in 0:1) {
	# Asymmetric is second
	rows1 <- which((asymm*1) == j) 
	for (i in 1:8) {
		rows2<- which (states[, i] == 1)
		rows <- intersect(rows1, rows2)
		if (j == 0 ) {
			symmMat[i, 2] <- mean(correct[rows])
		} else {
			asymmMat[i, 2] <- mean(correct[rows])
		}
	}
}
# by node


#  1 variable, 8 categories
findState <- function(row) {which(row != 0) }
truth <- apply(states, 1, findState) - 1
est <- read.csv("../8_cat/estimate/est2.empirical_est.labels_cat.csv")

colZero <- 2
colLast <- colZero + 7 

maxRow <- apply(est[,colZero:colLast], 1, max)
pointEst <- rep(NA, reps)
for (i in 1:reps) {
	pointEst[i] <- which(est[i, colZero:colLast] ==  maxRow[i]) - 1 
}

correct <- pointEst == truth 
#correct <- pointEst == truth[,2] 
print("correct 8")
mean(correct)

print("by states")
# These are comparing different things- individual nodes correct vs all nodes on a tree
for (j in 0:1) {
	# Asymmetric is second
	rows1 <- which((asymm*1) == j) 

	for (i in 0:7) {
		rows2<- which (truth  == i)
		rows <- intersect(rows1, rows2)
		if (j == 0 ) {
			symmMat[i +1, 3] <- mean(correct[rows])
		} else {
			asymmMat[i+ 1, 3] <- mean(correct[rows])
		}
	}
}

colnames(symmMat) <- c("marginal", "single node", "joint")
colnames(asymmMat) <- c("marginal", "single node", "joint")
write.csv(symmMat, "sym_train2.csv", row.names = FALSE)
write.csv(asymmMat, "asym_train2.csv", row.names = FALSE)
