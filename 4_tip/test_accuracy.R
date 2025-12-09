# 1 cateogorical with 8 states
true8 <- read.csv("~/asr_phyddle/4_tip/8_cat/estimate/out.test_true.labels_cat.csv")
est  <-read.csv("~/asr_phyddle/4_tip/8_cat/estimate/out.test_est.labels_cat.csv")
order <-read.csv("~/asr_phyddle/4_tip/3_bin/order.csv")

# Find what individuals nodes are based on how the nodes are coded
true8$node_1 <- 0
true8[which(true8$asr_0 >3), 3] <- 1

true8$node_2 <- 0
true8[which(true8$asr_0  ==2), 4] <- 1
true8[which(true8$asr_0  ==3), 4] <- 1
true8[which(true8$asr_0  ==6), 4] <- 1
true8[which(true8$asr_0  ==7), 4] <- 1

true8$node_3 <- 0
true8[which(true8$asr_0  ==1), 5] <- 1
true8[which(true8$asr_0  ==3), 5] <- 1
true8[which(true8$asr_0  ==5), 5] <- 1
true8[which(true8$asr_0  ==7), 5] <- 1

# Find marginal point estimates, which state has > .5 probability
nodeEst <- matrix(NA, nrow = 2500, ncol = 3)
nodeEst[, 1] <- apply(est[, 6:9], 1, sum) > .5
nodeEst[, 2] <- apply(est[,c(5,6, 8,9)], 1, sum) > .5
nodeEst[, 3] <- apply(est[,c(3, 5, 7,9)], 1, sum) > .5
nodeEst <- nodeEst * 1.0

# Find marginal probabilities
est8 <- matrix(NA, nrow = 2500, ncol = 3)
est8[, 1] <- apply(est[, 6:9], 1, sum)
est8[, 2] <- apply(est[,c(5,6, 8,9)], 1, sum) 
est8[, 3] <- apply(est[,c(3, 5, 7,9)], 1, sum) 

# This is the overall accuracy
mean(true8[, 3:5] == nodeEst)
est8Cat <- rep(NA, 7500)
true8Cat <- rep(NA, 7500)

for(i in 1:2500) {
  rows <- ((3 *(i-1)) + 1): ((3 *(i-1)) + 3)
  for (j in 1:3) {
    new <- which (order[rows, 2] == j)
    true8Cat[(3 *(i-1)) + j] <- true8[i , 2 + new ]
    est8Cat[(3 *(i-1)) + j] <- nodeEst[i, new]
  }
  
}

# 3 binary variables
est3 <- read.csv("~/asr_phyddle/4_tip/3_bin/estimate/out.test_est.labels_cat.csv")
true3 <- read.csv("~/asr_phyddle/4_tip/3_bin/estimate/out.test_true.labels_cat.csv")

# Find point estimate (which state > .5 probability)
est3_bin <- matrix(NA, nrow = 2500, ncol = 3)
est3_bin[, 1] <- (est3[, 3] > .5) * 1
est3_bin[, 2] <- (est3[, 5] > .5) * 1
est3_bin[, 3] <- (est3[, 7] > .5) * 1

est3Cat <- rep(NA, 7500)
for(i in 1:2500) {
  rows <- ((3 *(i-1)) + 1): ((3 *(i-1)) + 3)
  for (j in 1:3) {
    new <- which (order[rows, 2] == j)
    est3Cat[(3 *(i-1)) + j] <- est3_bin[i, new]
  }

}
# Accuracy of 3 binary estimates
mean(est3_bin == true3[, 2:4])

# true3 <- read.csv("~/asr_phyddle/4_tip/3_bin/all.anc_state.csv", header = FALSE)
# # Accuracy of 3 binary estimates, just in a different order
# mean(est3Cat == true3[, 2])

# ML
bayes <- read.csv("~/asr_phyddle/4_tip/3_bin/bayes/all_Bayes.csv")
probZero <- which(is.na(bayes$anc_state_2))
bayes$anc_state_2[probZero] <- (bayes$anc_state_1 == 0)[probZero]
row_Ones <- which(bayes$anc_state_1 == 1)
row_Twos <- which(bayes$anc_state_2 == 1)
bayes$probOne <- NA
bayes$probOne[row_Ones] <- bayes$anc_state_1_pp[row_Ones] 
bayes$probOne[row_Twos] <- bayes$anc_state_2_pp[row_Twos] 
ml3 <- bayes[, c(1:6)]
true3_AS <- read.csv("~/asr_phyddle/4_tip/3_bin/all.anc_state.csv", header = FALSE)
#ml3 <- read.csv("~/asr_phyddle/4_tip/3_bin/all.asr.csv", header = FALSE)

# reformat_ml3 <- cbind(ml3$V4[seq(1, to = 7500, by = 3)], 
#                       ml3$V4[seq(2, to = 7500, by = 3)], 
#                       ml3$V4[seq(1, to = 7500, by = 3)])

#diff_est_ml <- abs(est8 - reformat_ml3)
# pdf("~/asr_writing/figs/4_tip/one_cat.pdf", width=6, height = 4)
# hist(as.vector(diff_est_ml), xlim = c(0, .9), ylim = c(0, 4000), xlab = "abs(Phyddle probability - EB probability)", main = "One categorical, 8 states")
# dev.off()
# mean(diff_est_ml)
# 
# reformat_est3 <- as.matrix(est3[, c(5, 7, 9)])
# diff_est_ml <- abs(reformat_est3 - reformat_ml3)
# pdf("~/asr_writing/figs/4_tip/three_cat.pdf", width=6, height = 4)
# hist(as.vector(diff_est_ml), xlim = c(0, .9), ylim = c(0, 4000), xlab = "abs(Phyddle probability - EB probability)", main = "3 categorical variable")
# dev.off()
# mean(diff_est_ml)


mean(true3_AS[,2] == ml3[,2])

# 1 state at a time
trueOne   <- read.csv("~/asr_phyddle/4_tip/one/truth.csv", header = FALSE)
estOne    <- read.csv("~/asr_phyddle/4_tip/one/estimate/out.empirical_est.labels_cat.csv")
estOneBin <- (estOne[, 3] == apply(estOne[,2:3], 1, max))* 1

# Accuracy of 1 state at a time
mean(trueOne == estOneBin)

#match to ml 
print("Ml to one node at at time")
mean(ml3[,2] == estOneBin)
print("Ml to one categorical with 8 states")
mean(ml3[,2] == est8Cat)
print("Ml to 3 binary variables ")
mean(ml3[,2] == est3Cat)

# Ml accuracy
print("Ml overall accuracy")

mean(ml3[,2] ==  true3_AS[,2])

# Match to ML by node
node1 <- seq(1, by = 3, to = 7500)

mean(ml3[node1,2] == estOneBin[node1])
mean(est8Cat[node1] == ml3[node1,2])
mean(ml3[node1,2] == est3Cat[node1])

node2 <- seq(2, by = 3, to = 7500)
mean(ml3[node2,2] == estOneBin[node2])
mean(est8Cat[node2] == ml3[node2,2])
mean(ml3[node2,2] == est3Cat[node2])

node3 <- seq(3, by = 3, to = 7500)
mean(ml3[node3,2] == estOneBin[node3])
mean(est8Cat[node3] == ml3[node3,2])
mean(ml3[node3,2] == est3Cat[node3])

reformat_estOne <- cbind(estOne$asr_node_state_1[node1], 
      estOne$asr_node_state_1[node2], 
      estOne$asr_node_state_1[node3])
diff_est_ml <- abs(reformat_estOne - reformat_ml3)
# pdf("~/asr_writing/figs/4_tip/one_node.pdf", width=6, height = 4)
# hist(as.vector(diff_est_ml), xlim = c(0, .9), ylim = c(0, 4000), xlab = "abs(Phyddle probability - EB probability)", main = "One node at a time")
# dev.off()
mean(diff_est_ml)

# Matches so each other
print("One at a time to 8 cat percent match")
mean(estOneBin == est8Cat )
print("3 binary to 8 cat percent match")
mean( est3Cat == est8Cat )
print("One at a time to 3 binary percent match")
mean(estOneBin ==  est3Cat )

mean((true3[,2] == trueOne))
mean(true8Cat == trueOne)

print("8 cat to truth")
mean(true8Cat == est8Cat)
print("One at a time to truth")
mean(trueOne == estOneBin)
print("Accuracy of 3 binary estimates")
mean(est3_bin == true3[, 2:4])
