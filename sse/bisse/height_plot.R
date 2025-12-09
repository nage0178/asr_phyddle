library(caret)
library(wrapr)
library(ggplot2)
library("rhdf5")
library(cowplot)
library(ggpubr)

savedir <- "~/asr_writing/manuscript/figs/"

dirs <-  c("~/asr_phyddle/mk_binary/fix/50_unequal_q/", "~/asr_phyddle/sse/bisse_stat_q/")
numBins <- 11
numDirs <- 2
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
  
 # if (dir == dirs[2])  {
    bayes <- read.csv(paste(dir, "bayes/all_Bayes.csv", sep = ""))
    probZero <- which(is.na(bayes$anc_state_2))
    bayes$anc_state_2[probZero] <- (bayes$anc_state_1 == 0)[probZero]
    row_Ones <- which(bayes$anc_state_1 == 1)
    row_Twos <- which(bayes$anc_state_2 == 1)
    bayes$probOne <- NA
    bayes$probOne[row_Ones] <- bayes$anc_state_1_pp[row_Ones] 
    bayes$probOne[row_Twos] <- bayes$anc_state_2_pp[row_Twos] 
 # } else {
  #  ml <- read.csv(paste(dir, "ml_asr.csv", sep = ""), header = FALSE)
    
  #}

  
  nstates <- dim(truth)[2] - 1 
  
  # Binary estimates
  bin_est <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2])
  prob    <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  bin_phy <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  
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
  #if (dir == dirs[2]) {
    correct_est_ml <- (matrix(bayes$anc_state_1, nrow =2500, byrow = TRUE )== truth[, 2:dim(truth)[2]])
 # } else {
  #  correct_est_ml <- (ml == truth[, 2:dim(truth)[2]])
  #}

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
    print(length(index))
    mean_correct[i+1, 1] <- mean(times_by_est[index,1])
    mean_correct[i+1, 2] <- mean(times_by_est[index,2])
    
    mean_correct_ml[i+1, 1] <- mean(times_by_est[index,1])
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
}

method <- c(rep("phyddle", dim(correct_by_height)[1]), rep("Bayesian\ninference", dim(correct_by_height_ml)[1]))

correct_by_height_df <- as.data.frame(rbind(correct_by_height, correct_by_height_ml))
correct_by_height_df$n_tips <- as.factor(correct_by_height_df$n_tips )
correct_by_height_df$size_range <- as.factor(correct_by_height_df$size_range )
correct_by_height_df <- cbind(correct_by_height_df, method)

size_names <- c(
  `0` = "Markov",
  `1` = "BISSE")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(2)

#pdf(paste(savedir, "mk_height_fixVarSize.pdf", sep = ""), width =8, height = 3)
plot1 <- ggplot(correct_by_height_df, aes(avg_height, avg_correct, color = size_range, pch = method, linetype = method)) + #facet_wrap(~size_range)+ #, 
                                                                                                                                 #labeller = as_labeller(size_names)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean proportion of tree height\n ",
       y = "proportion correct", 
       color = "model") + 
  scale_color_manual(
    values = cols, # Assign specific colors
    labels = c("BiSSE", "Markov") # Assign custom labels
  )+ theme(legend.position = "none")

plotleg<- ggplot(correct_by_height_df, aes(avg_height, avg_correct, color = size_range, pch = method, linetype = method)) + #facet_wrap(~size_range)+ #, 
  #labeller = as_labeller(size_names)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean proportion of tree height\n ",
       y = "probability correct", 
       color = "model") + 
  scale_color_manual(
    values = cols, # Assign specific colors
    labels = c("BiSSE", "Markov") # Assign custom labels
  )+ theme(legend.position = "bottom")

legend1 <- cowplot::get_legend(plotleg)
leg1 <- as_ggplot(legend1)

probs_mat  <- matrix(NA, ncol = 2, nrow = 0)
colnames(probs_mat) <- c("phyddle", "Bayesian inference")

simple_estimator <- read.csv("~/asr_phyddle/sse/bisse_stat_q/prop_tips_clade_inference_correct.csv")
df <- as.data.frame(simple_estimator)
small <- correct_by_height_df[intersect(which(correct_by_height_df$method == "phyddle"), which(correct_by_height_df$size_range == "0")), ]
colnames(small)[2:3] <- colnames(df)[2:3]
df <- rbind(df[2:3], small[, 2:3])
df$method <- c(rep("naive estimator", 11), rep("phyddle", 11))
plot_supp <- ggplot(df, aes(avg_age, avg_correct,  pch = method, linetype = method)) + #facet_wrap(~size_range)+ #, 
  #labeller = as_labeller(size_names)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean proportion of tree height\n ",
       y = "proportion correct", 
       color = "model") 
pdf(paste(savedir, "bisse_naive_est.pdf", sep = ""), width =5, height = 4)
plot_supp
dev.off()

diff <- c()
for (dir in dirs) {
  # Reordered plots
  #if (dir == dirs[2]) {
    bayes <- read.csv(paste(dir, "bayes/all_Bayes.csv", sep = ""))
    probZero <- which(is.na(bayes$anc_state_2))
    bayes$anc_state_2[probZero] <- (bayes$anc_state_1 == 0)[probZero]
    row_Ones <- which(bayes$anc_state_1 == 1)
    row_Twos <- which(bayes$anc_state_2 == 1)
    bayes$probOne <- NA
    bayes$probOne[row_Ones] <- bayes$anc_state_1_pp[row_Ones]
    bayes$probOne[row_Twos] <- bayes$anc_state_2_pp[row_Twos]
    ml_prob <- matrix(bayes$probOne, nrow = 2500, byrow = TRUE)
  #} else {
  #  dir <- "~/asr_phyddle/mk_binary/fix/50_unequal/"
  #  ml_prob <- read.csv(paste(dir, "ML_prob.csv", sep = ""))
  #}
  est <- read.csv(paste(dir, "est.csv", sep = ""))
  
  if (  grepl("fix", dir, ignore.case = FALSE)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  
  #est_phyddle <- read.csv(paste(dir, "estimate/out.test_est.labels_cat.csv", sep = ""))
  estOne <- est[, seq(3, dim(est)[2], by = 2)]
  diff1 <- c(as.matrix(estOne-ml_prob))
  diff <- c(diff, diff1)
  probs <- cbind(c(as.matrix(estOne)), c(as.matrix(ml_prob)))
  probs_mat <- rbind(probs_mat, probs)
}
numRow <- dim(probs_mat)[1]
set.seed(1)
#reorder <- sample(1:numRow)[1:2500] 
nodeFromTree <- sample(1:49, 5000, replace = TRUE)
startIndex <- 0:4999 * 49 
reorder <- sample(nodeFromTree + startIndex)

probs_df <- as.data.frame(probs_mat)
probs_df$model <- c(rep("Markov", 2500 * 49), rep("BiSSE", 2500 * 49))
ggplot(probs_df[reorder, ], aes(x =`Bayesian inference`, y = phyddle, )) + geom_point()

print("correlations")
cor(probs_df$phyddle[1:(2500*49)], probs_df$`Bayesian inference`[1:(2500*49)])
cor(probs_df$phyddle[(2500*49 + 1):nrow(probs_df)], probs_df$`Bayesian inference`[(2500*49 + 1):nrow(probs_df)])


plot2 <- ggplot(probs_df[reorder[1:500],], aes(x =`Bayesian inference`, y = phyddle, color = model,)) + 
  geom_point(size = .1) + theme_classic()+ theme(legend.position = "none") + 
  labs(x = "Bayesian inference\n")

probs_df$diff <- probs_df$phyddle - probs_df$`Bayesian inference`

plot3 <- 
  ggplot(probs_df, aes(x = diff, fill = model)) + 
  geom_histogram(data=subset(probs_df, model == "Markov"), aes(y = after_stat(count / sum(count))) ,   bins = 100, alpha = .7) + 
  geom_histogram(data=subset(probs_df, model == "BiSSE"), aes(y = after_stat(count / sum(count))) ,   bins = 100, alpha = .7) + 
  theme_classic() + labs(x = "phyddle probability -\nBayesian inference probability", y = "proportion of inferences") +
  theme(legend.position = "none")
model <- c(rep("Markov", 2500 * 49), rep("BiSSE", 2500 * 49))
df <- data.frame (model, diff)
# plot2 <- ggplot(df, aes(x = diff, color = model, fill = model)) + geom_histogram(position = "identity",  alpha = 0.5,) + 
#   labs(x = "phyddle probability -\n maximum likelihood probability",
#        color = "model") + theme_classic() + theme(legend.position = "none")


mean(abs(probs_df[1:(dim(probs_df)[1]/2), 4]) > .5)
mean(abs(probs_df[1:(dim(probs_df)[1]/2), 4]) > .2585)
mk_int <- .2585
mean(abs(probs_df[((dim(probs_df)[1]/2)+1):dim(probs_df)[1], 4]) > .5)
mean(abs(probs_df[((dim(probs_df)[1]/2)+1):dim(probs_df)[1], 4]) > .386)
bisse_int <- .386

plot2_lines <- plot2 + geom_abline(intercept = bisse_int , slope = 1, color = cols[1], linetype = 2) + geom_abline(intercept = -bisse_int , slope = 1, color = cols[1], linetype = 2) + 
  geom_abline(intercept = mk_int , slope = 1, color = cols[2], linetype = 2) + geom_abline(intercept = -mk_int , slope = 1, color = cols[2], linetype = 2) 

pdf(paste(savedir, "bisse.pdf", sep = ""), width =8, height = 3)
plotSet <- ggarrange(plot1, plot3, plot2_lines,  nrow =1 , ncol = 3, labels= c("a", "b", "c"), widths = c(.33 * 8, .33 * 8 , .33 * 8 ))
ggarrange(plotSet, leg1, nrow = 2, ncol = 1, widths = 8, heights = c(3 * .9, 3 * .1))
dev.off()

png(paste(savedir, "bisse_a.png", sep = ""), width =5, height = 3, units ="in", res = 500)
ggplot(correct_by_height_df, aes(avg_height, avg_correct, color = size_range, pch = method, linetype = method)) + #facet_wrap(~size_range)+ #, 
  #labeller = as_labeller(size_names)) + 
  geom_point() + theme_classic() + geom_line() + 
  labs(x = "mean proportion of tree height\n ",
       y = "proportion correct", 
       color = "model") + 
  scale_color_manual(
    values = cols, # Assign specific colors
    labels = c("BiSSE", "Markov") # Assign custom labels
  )
dev.off()

png(paste(savedir, "bisse_c.png", sep = ""), width =5, height = 3, units ="in", res = 500)
ggplot(probs_df[reorder[1:500],], aes(x =`Bayesian inference`, y = phyddle, color = model,)) + 
  geom_point(size = .1) + theme_classic()+
  labs(x = "Bayesian inference\n") + coord_fixed(ratio = 1)
dev.off()

png(paste(savedir, "bisse_c2.png", sep = ""), width =5, height = 3, units ="in", res = 500)
ggplot(probs_df[reorder[1:500],], aes(x =`Bayesian inference`, y = phyddle, color = model,)) + 
  geom_point(size = .1) + theme_classic()+coord_fixed(ratio = 1) +
  labs(x = "Bayesian inference\n") + geom_abline(intercept = bisse_int , slope = 1, color = cols[1], linetype = 2) + geom_abline(intercept = -bisse_int , slope = 1, color = cols[1], linetype = 2) + 
  geom_abline(intercept = mk_int , slope = 1, color = cols[2], linetype = 2) + geom_abline(intercept = -mk_int , slope = 1, color = cols[2], linetype = 2) 
dev.off()

png(paste(savedir, "bisse_b.png", sep = ""), width =5, height = 3, units ="in", res = 500)
ggplot(probs_df, aes(x = diff, fill = model)) + 
  geom_histogram(data=subset(probs_df, model == "Markov"), aes(y = after_stat(count / sum(count))) ,   bins = 100, alpha = .7) + 
  geom_histogram(data=subset(probs_df, model == "BiSSE"), aes(y = after_stat(count / sum(count))) ,   bins = 100, alpha = .7) + 
  theme_classic() + labs(x = "phyddle probability -\nBayesian inference probability", y = "proportion of inferences") 
dev.off()

truth <- read.csv(paste(dirs[1], "truth.csv", sep = ""))
# hist(apply(truth[, 2:dim(truth)[2]], 1, mean), xlab = "frequency of state 1", main="Binary Markov model")
# truth <- read.csv(paste(dirs[2], "truth.csv", sep = ""))
# hist(apply(truth[, 2:dim(truth)[2]], 1, mean), xlab = "frequency of state 1", main="BISSE")
# table(apply(truth[, 2:dim(truth)[2]], 1, mean))
# truth <- read.csv("~/asr_phyddle/mk_binary/fix/50_unequal/truth.csv")
# hist(apply(truth[, 2:dim(truth)[2]], 1, mean), xlab = "frequency of state 1", main="Mk unequal")
zero_one <- union (which(probs_df$phyddle > 1 - 10^-5), which(probs_df$phyddle <  10^-5))
cor(probs_df$phyddle[-zero_one], probs_df$`Bayesian inference`[-zero_one])
ggplot(probs_df[-zero_one], aes(x =`Bayesian inference`, y = phyddle)) + geom_point(size = .01)
