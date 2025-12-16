library(ggplot2)
dir <- "~/asr_phyddle/4_tip/3_bin/"
truth <- read.csv(paste(dir, "truth.csv", sep = ""))
est <- read.csv(paste(dir, "est.csv", sep = ""))
ml <- read.csv(paste(dir, "ml_asr.csv", sep = ""), header = FALSE)
ml_prob <- as.matrix(read.csv(paste(dir, "ml_asr_prob.csv", sep = ""), header = FALSE))
est_1 <- cbind(est$asr_1_1, est$asr_1_2, est$asr_1_3)

df <- cbind(c(ml_prob), c(est_1))
colnames(df) <- c("maximum likelihood", "phyddle")
df <- as.data.frame(df)
ggplot(df, aes(x =`maximum likelihood`, y = phyddle )) + geom_point()

dir <- "~/asr_phyddle/4_tip/one/"

est_one <- read.csv(paste(dir, "estimate/out.empirical_est.labels_cat.csv", sep = ""))
df_one <- cbind(c(t(ml_prob)), est_one$asr_node_state_1)
colnames(df_one) <- c("maximum likelihood", "phyddle")
df_one <- as.data.frame(df_one)
ggplot(df_one, aes(x =`maximum likelihood`, y = phyddle )) + geom_point()
cor(df_one$`maximum likelihood`, df_one$phyddle)

df_both <- rbind(df, df_one)
df_both$method <- c(rep("marginal", 7500), rep("single node", 7500))
reorder <- sample(1:15000) 
ggplot(df_both[reorder, ], aes(x =`maximum likelihood`, y = phyddle, color = method )) + geom_point() + theme_classic()


bayes <- read.csv("~/asr_phyddle/4_tip/3_bin/bayes/all_Bayes.csv")
# If probability is na, that means the probability is zero
probZero <- which(is.na(bayes$anc_state_2))
bayes$anc_state_2[probZero] <- (bayes$anc_state_1 == 0)[probZero]
# Find which rows prefer state 1 vs 2
row_Ones <- which(bayes$anc_state_1 == 1)
row_Twos <- which(bayes$anc_state_2 == 1)
# Probability of most probable state
bayes$probOne <- NA
bayes$probOne[row_Ones] <- bayes$anc_state_1_pp[row_Ones] 
bayes$probOne[row_Twos] <- bayes$anc_state_2_pp[row_Twos] 
#
bayes$phyddle <- est_one$asr_node_state_1
bayes_both <- rbind(bayes, bayes, bayes)
# Marginal estimates
bayes_both$phyddle[1:dim(bayes)[1]] <- c(t(est_1))
bayes_both$method <- c(rep("marginal", 7500), rep("single node", 7500), rep("joint", 7500))
colnames(bayes_both)[6] <- "Bayesian"
#ggplot(bayes, aes(x=probOne, y = phyddle)) + geom_point() + theme_classic() + labs(x= "posterior probability", y = "phyddle probability")

dir <- "~/asr_phyddle/4_tip/8_cat/"
true8 <- read.csv(paste(dir, "estimate/out.test_true.labels_cat.csv", sep = ""))
est  <-read.csv(paste(dir, "estimate/out.test_est.labels_cat.csv", sep = ""))
order_mat <-read.csv(paste(dir, "../3_bin/order.csv", sep = ""))


# Find marginal point estimates, which state has > .5 probability
nodeEst <- matrix(NA, nrow = 2500, ncol = 3)
nodeEst[, 1] <- apply(est[, 6:9], 1, sum) #> .5
nodeEst[, 2] <- apply(est[,c(5,6, 8,9)], 1, sum)# > .5
nodeEst[, 3] <- apply(est[,c(3, 5, 7,9)], 1, sum)# > .5
#nodeEst <- nodeEst * 1.0

reorderEst <- matrix(NA, nrow=2500, ncol = 3)
for (i in 1:2500)  {
  rows <- (i * 3 - 2) : (i * 3)
  
  small <- order_mat[rows,]
  reorder <- (small[order(small$original), ])
  reorderEst[i, ] <- nodeEst[i, reorder$new + 1]
  #nodeEst[i, ] <- nodeEst[i, ]
}
bayes_both$phyddle[(dim(bayes)[1]*2 + 1):(dim(bayes)[1]*3) ] <- c(t(est_1))


set.seed(1)
node <- sample(1:3, 7500, replace = TRUE)
startIndex <- 0:7499* 3 
reorder <- sample(node + startIndex)
#reorder <- sample(1:(7500 * 3))[1:5000] 
pdf("~/asr_writing/manuscript/figs/4_tip_prob.pdf", width = 4, height = 3)
ggplot(bayes_both[reorder, ], aes(x =Bayesian, y = phyddle, color = method )) + geom_point(size = .1) + theme_classic()
dev.off()

png("~/asr_writing/manuscript/figs/4_tip_prob.png", width = 4, height = 3, units = "in", res = 500)
ggplot(bayes_both[reorder, ], aes(x =Bayesian, y = phyddle, color = method )) + geom_point(size = .1) + theme_classic()
dev.off()


# Match to Bayesian 
both1 <- which((bayes$phyddle > .5 ) * (bayes$probOne > .5)  ==1)
both0 <- which((bayes$phyddle < .5 ) * (bayes$probOne  < .5)  ==1)
(length(both1) + length(both0)) / 7500. 
