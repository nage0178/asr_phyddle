library(caret)
library(wrapr)
library(ggplot2)
setwd("~/asr_phyddle/mk_binary/")
wkdirs <- c("fix/", "var/")
savedir <- "~/asr_writing/manuscript/figs/"
dirs <- c("50/", "100/", "200/")
dirs <- c(paste(wkdirs[1], dirs, sep = ""), paste(wkdirs[2], dirs, sep = ""))

numBins <- 11
numDirs <- 6
correct_by_height <- matrix(NA, ncol = 4, nrow =numBins * numDirs )
colnames(correct_by_height) <- c("n_tips",  "avg_height", "avg_correct", "size_range")

correct_by_height_ml <- matrix(NA, ncol = 4, nrow =numBins * numDirs )
colnames(correct_by_height_ml) <- c("n_tips",  "avg_height", "avg_correct", "size_range")

probs_correct <- matrix(NA, ncol = 4, nrow =numBins * numDirs )
colnames(probs_correct) <- c("average_estimate", "probability_correct", "n_tips", "size_range")

j <- 1
for (dir in dirs) {
  if (  grepl("fix", dir, ignore.case = FALSE)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  truth <- read.csv(paste(dir, "truth.csv", sep = ""))
  est <- read.csv(paste(dir, "est.csv", sep = ""))
  ages <- read.csv(paste(dir, "prop_ages.csv", sep = ""))
  est_phyddle <- read.csv(paste(dir, "estimate/out.test_est.labels_cat.csv", sep = ""))
  true_phyddle <- read.csv(paste(dir, "estimate/out.test_true.labels_cat.csv", sep = ""))
  
  bayes <- read.csv(paste(dir, "bayes/all_Bayes.csv", sep = ""))
  conv  <- read.table(paste(dir,"bayes/notConverged", sep = ""), header = TRUE, sep = ",")
  idx_rm <- conv[, 2]
  #probZero <- which(is.na(bayes$anc_state_2))
  #bayes$anc_state_2[probZero] <- (bayes$anc_state_1 == 0)[probZero]
  #row_Ones <- which(bayes$anc_state_1 == 1)
  #row_Twos <- which(bayes$anc_state_2 == 1)
  #bayes$probOne <- NA
  #bayes$probOne[row_Ones] <- bayes$anc_state_1_pp[row_Ones]
  #bayes$probOne[row_Twos] <- bayes$anc_state_2_pp[row_Twos]
  #ml <- matrix(bayes$probOne, nrow = 2500, byrow = TRUE)
  ml <- matrix(bayes$anc_state_1, ncol = dim(truth)[2] -1, byrow = TRUE)
  ml[idx_rm, ] <- NA
  
  if (dim(est_phyddle)[1] > 2500) {
    est_phyddle <- est_phyddle[1:2500, ]
    true_phyddle <- true_phyddle[1:2500, ]
  }
  
  nstates <- dim(truth)[2] - 1 
  
  #results <-matrix(NA, nrow = nstates, ncol = 7)
  
  # Binary estimates
  bin_est <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2])
  prob    <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  bin_phy <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  
  #colnames(results) <- c("DL", "ML", "DL_ML", "DL_low", "DL_high", "ML_low", "ML_high")
  for (i in c(1:(dim(est_phyddle)[2]/2))) {
    colZero <- 2+(i-1)*2
    colOne <- 3+(i-1)*2

    # Most likely state is 1 
    infer <- (est[,colOne] == apply(cbind(est[,colZero], est[,colOne] ), 1, max)) * 1
    bin_est[, i+1] <- infer
    
    if (colOne <= dim(est_phyddle)[2]) {
      infer <- (est_phyddle[,colOne] == apply(est_phyddle[,colZero:colOne], 1, max)) * 1
      prob[, i] <- est_phyddle[,colOne]
    } else {
      infer <- 1 - est_phyddle[,colZero]
      prob[, i] <-  1- est_phyddle[,colZero]
    }
    
    bin_phy[,i] <- 1 == factor(true_phyddle[,i + 1])

  }
  
  correct_prop <- cbind(c(prob), c(bin_phy), c(floor(prob* 10)))
  plot_prop <- matrix(NA, nrow = 11, ncol = 2)
  
  for (i in 0:10) {
    index <- which(correct_prop[,3] == i)
    #  print(length(index))
    plot_prop[i+1, 1] <- mean(correct_prop[index,1])
    plot_prop[i+1, 2] <- mean(correct_prop[index,2])
  }
  
  indeces <- ((j-1) *numBins +1): (j*numBins)
  probs_correct[indeces, 3] <- nstates + 1
  probs_correct[indeces, 1:2] <- plot_prop[, 1:2]
  probs_correct[indeces, 4] <- fix
  
  # This is the proportion correct by tree height
  correct_est <- (bin_est == truth)[,2:dim(truth)[2]]
  correct_est_ml <- (ml == truth[, 2:dim(truth)[2]])
  correct_est_ml[idx_rm, ] <- NA
 
  v_ages <- c(as.matrix(ages[,2:dim(truth)[2]]))
  times_by_est <- cbind(v_ages, c(correct_est),  c(correct_est_ml))
  colnames(times_by_est) <- c("age", "phyddle_correct", "ml_correct")
  # bc trees are of variable size
  if (length(which(is.na(times_by_est[,2]))) != 0 ) {
    times_by_est <- times_by_est[-which(is.na(times_by_est[,2])),]
  }
  
  times_by_est <- cbind(times_by_est, floor(times_by_est[,1] * 10 ))
  colnames(times_by_est)[4] <- "bin"
  mean_correct <- matrix(0, nrow = 11, ncol = 2)
  colnames(mean_correct) <- c("avg_age", "avg_phyddle_correct")
  
  mean_correct_ml <- matrix(0, nrow = 11, ncol = 2)
  colnames(mean_correct_ml) <- c("avg_age", "avg_ml_correct")
  
  for (i in 0:10) {
    index <- which(times_by_est[,4] == i)
    
    mean_correct[i+1, 1] <- mean(times_by_est[index,1])
    mean_correct[i+1, 2] <- mean(times_by_est[index,2])
    
    mean_correct_ml[i+1, 1] <- mean(times_by_est[index,1], na.rm = TRUE)
    mean_correct_ml[i+1, 2] <- mean(times_by_est[index,3], na.rm = TRUE)
  }
  
  correct_by_height[((j-1) *numBins +1): (j*numBins), 1] <- nstates + 1
  correct_by_height[((j-1) *numBins +1): (j*numBins), 2] <- mean_correct[,1]
  correct_by_height[((j-1) *numBins +1): (j*numBins), 3] <- mean_correct[,2]
  correct_by_height[((j-1) *numBins +1): (j*numBins), 4] <- fix
  
  correct_by_height_ml[((j-1) *numBins +1): (j*numBins), 1] <- nstates + 1
  correct_by_height_ml[((j-1) *numBins +1): (j*numBins), 2] <- mean_correct_ml[,1]
  correct_by_height_ml[((j-1) *numBins +1): (j*numBins), 3] <- mean_correct_ml[,2]
  correct_by_height_ml[((j-1) *numBins +1): (j*numBins), 4] <- fix

  j <- j + 1
  print(paste("mean", dir))
  print(mean(bin_est[,2:dim(bin_est)[2]] == truth[,2:dim(bin_est)[2]], na.rm = TRUE))
}

method <- c(rep("phyddle", dim(correct_by_height)[1]), rep("Bayesian\ninference", dim(correct_by_height_ml)[1]))

correct_by_height_df <- as.data.frame(rbind(correct_by_height, correct_by_height_ml))
correct_by_height_df$n_tips <- as.factor(correct_by_height_df$n_tips )
correct_by_height_df$size_range <- as.factor(correct_by_height_df$size_range )
correct_by_height_df <- cbind(correct_by_height_df, method)

size_names <- c(
  `0` = "variable",
  `1` = "fixed")

pdf(paste(savedir, "mk_height_fixVarSize.pdf", sep = ""), width =8, height = 3)
both <- ggplot(correct_by_height_df, aes(avg_height, avg_correct, color = n_tips, pch = method, linetype = method)) + facet_wrap(~size_range, labeller = as_labeller(size_names)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean proportion of tree height",
       y = "probability inferred state is correct", 
       color = "tree size") 
print(both)
dev.off()

correct_by_height_fix <- correct_by_height_df[-which(correct_by_height_df$size_range == 0  ), ]
panelA <- ggplot(correct_by_height_fix, aes(avg_height, avg_correct, color = n_tips, pch = method, linetype = method)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "Mean proportion of tree height",
       y = "proportion correct", 
       color = "tree size")  +coord_cartesian(xlim = c(0, 1.01), ylim = c(.5, 1), expand = FALSE)  #ylim(0, 1)
pdf(paste(savedir, "mk_height_fixSize.pdf", sep = ""), width =5, height = 3)
print(panelA )
dev.off()

png(paste(savedir, "mk_height_fixSize.png", sep = ""), width =5, height = 3, units = "in", res = 500)
print(panelA )
dev.off()

correct_by_height_var <- correct_by_height_df[which(correct_by_height_df$size_range == 0  ), ]
pdf(paste(savedir, "mk_height_varSize.pdf", sep = ""), width =5, height = 3)
accuracy_height <- ggplot(correct_by_height_var, aes(avg_height, avg_correct, color = n_tips, pch = method, linetype = method)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "Mean proportion of tree height",
       y = "proportion correct", 
       color = "tree size") 
print(accuracy_height)
dev.off()

plot <- ggplot(correct_by_height_df, aes(avg_height, avg_correct, color = n_tips, pch = size_range, linetype = size_range)) + facet_wrap(~method) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean proportion of tree height",
       y = "probability inferred state is correct", 
       color = "tree size") 

probs_correct_df <- as.data.frame(probs_correct)
probs_correct_df$n_tips <- as.factor(probs_correct_df$n_tips)
probs_correct_df <- probs_correct_df[1:(dim(probs_correct_df)[1]/2), ]

pdf("~/asr_writing/manuscript/figs/prob_correct_by_prob.pdf", width =4, height = 3)
probCorrect <- ggplot(probs_correct_df, aes(average_estimate, probability_correct, color = n_tips, pch = n_tips, linetype = n_tips)) +
 # facet_wrap(~size_range, labeller = as_labeller(size_names)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean phyddle probability of state 1 ",
       y = "probability state 1 is correct", 
       color = "tree size", pch = "tree size", linetype = "tree size") 
print(probCorrect)
dev.off()

png("~/asr_writing/manuscript/figs/prob_correct_by_prob.png", width =4, height = 3, units = "in", res = 500)
print(probCorrect)
dev.off()
