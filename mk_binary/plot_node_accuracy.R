library(caret)
library(wrapr)
library(ggplot2)

wkdirs <- c("~/asr_phyddle/mk_binary/fix/", "~/asr_phyddle/mk_binary/var/")
savedir <- "~/asr_writing/manuscript/figs/"
dirs <- c("50/", "100/", "500/")
dirs <- c(paste(wkdirs[1], dirs, sep = ""), paste(wkdirs[2], dirs, sep = ""))

j <- 1
k <- 1
prob_by_node <- matrix(NA, ncol = 4, nrow = (49 + 99 + 499)*2)
colnames(prob_by_node) <- c("node_number", "probability_correct", "n_tips", "size_range")

for (dir in dirs) {
  truth <- read.csv(paste(dir, "truth.csv", sep = ""))
  est <- read.csv(paste(dir, "est.csv", sep = ""))
  ml <- read.csv(paste(dir, "ml_asr.csv", sep = ""), header = FALSE)
  
  if (  grepl("fix", dir, ignore.case = FALSE)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  
  # Reordered plots
  ml_reorder <- read.csv(paste(dir, "ml_est_reorder.csv", sep = ""))
  est_phyddle <- read.csv(paste(dir, "estimate/out.test_est.labels_cat.csv", sep = ""))
  true_phyddle <- read.csv(paste(dir, "estimate/out.test_true.labels_cat.csv", sep = ""))
  
  nstates <- dim(truth)[2] - 1 
  
  results <-matrix(0, nrow = nstates, ncol = 3)
  results_phyd <-matrix(0, nrow = nstates, ncol = 3)
  
  # Binary estimates
  bin_est <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2])
  prob <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  bin_phy <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  
  colnames(results) <- c("DL", "ML", "DL_ML")
  colnames(results_phyd) <- c("DL", "ML", "DL_ML")
  
  for (i in c(1:nstates)) {
    colZero <- 2+(i-1)*2
    colOne <- 3+(i-1)*2
    
    infer <- (est[,colOne] == apply(est[,colZero:colOne], 1, max)) * 1
    bin_est[, i+1] <- infer
    DL <-confusionMatrix( factor(infer), factor(truth[,i + 1]))
    ML <- confusionMatrix(factor(ml[,i]), factor(truth[,i +1] ))
  
    results[i,] <- c(DL$overall[1], ML$overall[1], mean(ml[,i] == infer, na.rm = TRUE))

    nodeExists <- which(!is.na(est[,colOne]))
    
    # For phyddle ordering, you are going to have to figure out what to fill in with NA
    infer <- (est_phyddle[,colOne] == apply(est_phyddle[,colZero:colOne], 1, max)) * 1
    prob[, i] <- est_phyddle[,colOne]
    bin_phy[,i] <- 1 == factor(true_phyddle[,i + 1])
    DL <- confusionMatrix( factor(infer[nodeExists]), factor(true_phyddle[nodeExists,i + 1]))
    ML <- confusionMatrix(factor(ml_reorder[,i+1]), factor(true_phyddle[,i +1] ))
    
    results_phyd[i,] <- c(DL$overall[1], ML$overall[1], mean(ml_reorder[,i+1] == infer) ) 

    prob_by_node[k, ] <- c(i, DL$overall[1], nstates + 1, fix)
    k <- k + 1
  }
  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(2)

prob_by_node_df <- as.data.frame(prob_by_node)
prob_by_node_df$n_tips <- as.factor(prob_by_node_df$n_tips)
prob_by_node_df$size_range <- as.factor(prob_by_node_df$size_range)

size_names <- c(
  `0` = "variable",
  `1` = "fixed")

pdf("~/asr_writing/manuscript/figs/accuracy_by_node.pdf", width =9, height = 3)
ggplot(prob_by_node_df, aes(node_number, probability_correct, color = size_range)) + geom_point(size = 0.5) + 
  facet_wrap(~n_tips, scales = "free_x", nrow = 1) + theme_classic()+  scale_color_manual(labels = c( "variable", "fixed"), values = cols) +
  labs(x = "Phyddle node number",
       y = "Accuracy", color = "tree size") 
dev.off()
