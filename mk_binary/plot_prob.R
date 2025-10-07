library(caret)
library(wrapr)
library("rhdf5")

wkdir <- "~/phyddle/workspace/fixSize/"
dirs <- c("32_tip/", "100_tip/", "500_tip/")
dirs <- paste(wkdir, dirs, sep = "")


j <- 1
numBins <- 11
probs_correct <- matrix(NA, ncol = 3, nrow =numBins * 3)
colnames(probs_correct) <- c("average_estimate", "probability_correct", "n_tips")

for (dir in dirs) {
  truth <- read.csv(paste(dir, "truth.csv", sep = ""))
  est <- read.csv(paste(dir, "est.csv", sep = ""))
  ml <- read.csv(paste(dir, "ml_asr.csv", sep = ""), header = FALSE)
  rates <- read.csv(paste(dir, "rates.csv", sep = ""))
  ages <- read.csv(paste(dir, "prop_ages.csv", sep = ""))
  
  fileTest <- H5Fopen(paste(dir, "format/out.test.hdf5", sep = ""))
  real_age <- fileTest$aux_data[2,]
  H5Fclose(fileTest)

  # Reordered plots
  ml_reorder <- read.csv(paste(dir, "ml_est_reorder.csv", sep = ""))
  est_phyddle <- read.csv(paste(dir, "estimate/out.test_est.labels_cat.csv", sep = ""))
  true_phyddle <- read.csv(paste(dir, "estimate/out.test_true.labels_cat.csv", sep = ""))
  
  # This is because I am currently reporting one "other" discrete character
  true_phyddle <- true_phyddle[,-2]
  est_phyddle <- est_phyddle[, -c(2:3)]
  
  nstates <- dim(truth)[2] - 1 
  low <- which(rates[,2] < 0.0277)
  high <- which(rates[,2] > 0.0277)
  
  old <- which (real_age  > 1.68 )
  young <- which (real_age <= 1.68)
  
  results <-matrix(0, nrow = nstates, ncol = 7)
  results_phyd <-matrix(0, nrow = nstates, ncol = 7)
  
  # Binary estimates
  bin_est <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2])
  prob <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  bin_phy <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  
 # colnames(results) <- c("DL", "ML", "DL_ML", "DL_low", "DL_high", "ML_low", "ML_high")
 # colnames(results_phyd) <- c("DL", "ML", "DL_ML", "DL_young", "DL_old", "ML_young", "ML_old")
  
  for (i in c(1:nstates)) {
    colZero <- 2+(i-1)*2
    colOne <- 3+(i-1)*2
   # print(colnames(est)[2+(i-1)*2])
   # print(colnames(est)[3+(i-1)*2])
    
    infer <- (est[,colOne] == apply(est[,colZero:colOne], 1, max)) * 1
    bin_est[, i+1] <- infer
    DL <-confusionMatrix( factor(infer), factor(truth[,i + 1]))
    ML <- confusionMatrix(factor(ml[,i]), factor(truth[,i +1] ))
    
    DL_low <-confusionMatrix( factor(infer[low]), factor(truth[low,i + 1]))
    DL_high <-confusionMatrix( factor(infer[high]), factor(truth[high,i + 1]))
    
    ML_low <- confusionMatrix(factor(ml[low,i]), factor(truth[low,i +1] ))
    ML_high <- confusionMatrix(factor(ml[high,i]), factor(truth[high,i +1] ))
    
  
    results[i,] <- c(DL$overall[1], ML$overall[1], mean(ml[,i] == infer, na.rm = TRUE), DL_low$overall[1], DL_high$overall[1], ML_low$overall[1], ML_high$overall[1])
    
    #print(DL$overall[1])
    #print(ML$overall[1])
    
    # For phyddle ordering, you are going to have to figure out what to fill in with NA
    infer <- (est_phyddle[,colOne] == apply(est_phyddle[,colZero:colOne], 1, max)) * 1
    prob[, i] <- est_phyddle[,colOne]
    bin_phy[,i] <- 1 == factor(true_phyddle[,i + 1])
    DL <-confusionMatrix( factor(infer), factor(true_phyddle[,i + 1]))
    ML <- confusionMatrix(factor(ml_reorder[,i+1]), factor(true_phyddle[,i +1] ))
    
    DL_old <-confusionMatrix( factor(infer[old]), factor(true_phyddle[old,i + 1]))
    DL_young <-confusionMatrix( factor(infer[young]), factor(true_phyddle[young,i + 1]))
    
    ML_old <- confusionMatrix(factor(ml_reorder[old,i+1]), factor(true_phyddle[old,i +1] ))
    ML_young <- confusionMatrix(factor(ml_reorder[young,i+1]), factor(true_phyddle[young,i +1] ))
    
    results_phyd[i,] <- c(DL$overall[1], ML$overall[1], mean(ml_reorder[,i+1] == infer), 
                          DL_young$overall[1], DL_old$overall[1], ML_young$overall[1], ML_old$overall[1])

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
  
 # pdf(paste("~/asr_writing/figs/fixSize/prop_correct_by_Prob", nstates + 1, ".pdf", sep = ""), width =8.67,  height = 8.67)
  plot(plot_prop[,1], plot_prop[,2], xlab = "average estimate of state 1 from phyddle", ylab = "probability 1 is true state", main = nstates+ 1)
  abline(a=0, b = 1)
 # dev.off()
  
  j <- j + 1
}
probs_correct_df <- as.data.frame(probs_correct)
probs_correct_df$n_tips <- as.factor(probs_correct_df$n_tips)
pdf("~/asr_writing/manuscript/figs/tmp_prob_correct_by_prob.pdf", width =4, height = 3)
ggplot(probs_correct_df, aes(average_estimate, probability_correct, color = n_tips, pch = n_tips, linetype = n_tips)) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean probability assigned to state 1 by phyddle",
       y = "probability state 1 is correct", 
       color = "tree size", pch = "tree size", linetype = "tree size") 
dev.off()
