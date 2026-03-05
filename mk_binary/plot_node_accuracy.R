library(caret)
library(wrapr)
library(ggplot2)

wkdirs <- c("~/asr_phyddle/mk_binary/fix/", "~/asr_phyddle/mk_binary/var/")
savedir <- "~/asr_writing/manuscript/figs/"
dirs <- c("50/", "100/", "200/")
dirs <- c(paste(wkdirs[1], dirs, sep = ""), paste(wkdirs[2], dirs, sep = ""))

j <- 1
k <- 1
prob_by_node <- matrix(NA, ncol = 4, nrow = (49 + 99 + 199)*2)
colnames(prob_by_node) <- c("node_number", "probability_correct", "n_tips", "size_range")

for (dir in dirs) {
  truth <- read.csv(paste(dir, "truth.csv", sep = ""))
  est <- read.csv(paste(dir, "est.csv", sep = ""))

  if (  grepl("fix", dir, ignore.case = FALSE)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  
  # Reordered plots
  est_phyddle <- read.csv(paste(dir, "estimate/out.test_est.labels_cat.csv", sep = ""))
  est_phyddle <- est_phyddle[1:2500, ]
  true_phyddle <- read.csv(paste(dir, "estimate/out.test_true.labels_cat.csv", sep = ""))
  true_phyddle <- true_phyddle[1:2500, ]
  
  nstates <- dim(truth)[2] - 1 
  
  # Binary estimates
  bin_est <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2])
  prob <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)
  bin_phy <- matrix(NA, nrow = dim(truth)[1], ncol = dim(truth)[2] -1)

  
  for (i in c(1:nstates)) {
    colZero <- 2+(i-1)*2
    colOne <- 3+(i-1)*2
    
    infer <- (est[,colOne] == apply(est[,colZero:colOne], 1, max)) * 1
    bin_est[, i+1] <- infer
    DL <-confusionMatrix( factor(infer), factor(truth[,i + 1]))

    nodeExists <- which(!is.na(est[,colOne]))
    
    # For phyddle ordering, you are going to have to figure out what to fill in with NA
    infer <- (est_phyddle[,colOne] == apply(est_phyddle[,colZero:colOne], 1, max)) * 1
    prob[, i] <- est_phyddle[,colOne]
    bin_phy[,i] <- 1 == factor(true_phyddle[,i + 1])
    DL <- confusionMatrix( factor(infer[nodeExists]), factor(true_phyddle[nodeExists,i + 1]))


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

panelC <- ggplot(prob_by_node_df, aes(node_number, probability_correct, color = size_range)) + geom_point(size = 0.5) + 
  facet_wrap(~n_tips, scales = "free_x", nrow = 1) + theme_classic()+  scale_color_manual(labels = c( "variable", "fixed"), values = cols) +
  labs(x = "phyddle node number",
       y = "proportion correct", color = "tree size") +coord_cartesian( ylim = c(0.5, 1), expand = FALSE)
pdf("~/asr_writing/manuscript/figs/accuracy_by_node.pdf", width =9, height = 3)
print(panelC)
dev.off()

png("~/asr_writing/manuscript/figs/accuracy_by_node.png", width =9, height = 3, units = "in", res = 500)
print(panelC)
dev.off()
