library(caret)
library(wrapr)
library(ggplot2)

wkdirs <- c("~/asr_phyddle/mk_binary/")
savedir <- "~/asr_writing/manuscript/figs/"
dirs <- c("50/", "100/", "200/")
dirs <- c(paste(wkdirs[1], dirs, sep = ""), paste(wkdirs[2], dirs, sep = ""))

true_file <- c("fix/50/estimate/out.test_true.labels_cat.csv", 
               "fix/50/estimate/out.test_true.labels_cat.csv",
               "fix/100/estimate/out.test_true.labels_cat.csv", 
               #"fix/100/estimate/out.test_true.labels_cat.csv", 
              # "fix/500/estimate/out.test_true.labels_cat.csv", 
               "fix/50/estimate/out.test_true.labels_cat.csv",
               "fix/100/estimate/out.test_true.labels_cat.csv", 
              
              "fix/50/estimate/out.test_true.labels_cat.csv",
              "fix/100/estimate/out.test_true.labels_cat.csv", 
              "fix/200/estimate/out.test_true.labels_cat.csv", "fix/200/estimate/out.test_true.labels_cat.csv"  )
# All of these need to be empirical
est_file <- c("var/50/estimate/50.empirical_est.labels_cat.csv", 
              "var/100/estimate/50.empirical_est.labels_cat.csv", 
              "var/100/estimate/100.empirical_est.labels_cat.csv", 
              #"var/500/estimate/100.empirical_est.labels_cat.csv",
              #"var/500/estimate/500.empirical_est.labels_cat.csv", 
              "fix/50/estimate/out.test_est.labels_cat.csv",
              "fix/100/estimate/out.test_est.labels_cat.csv", 
              
              "var/200/estimate/50.empirical_est.labels_cat.csv", 
              "var/200/estimate/100.empirical_est.labels_cat.csv", 
              "var/200/estimate/200.empirical_est.labels_cat.csv", "fix/200/estimate/out.test_est.labels_cat.csv" 
            )
tree_size <- c(50, 50, 100, #100, 500, 
               50, 100, 
               50, 100, 200, 200)
network_size <- c(50, 100, 100, #500, 500, 
                  50, 100, 
                  200, 200, 200, 200)
prob_by_node <- matrix(NA, ncol = 5, nrow = sum(tree_size - 1))
colnames(prob_by_node) <- c("node_number", "probability_correct", "n_tips", "network_size", "fixed")

j <- 1
k <- 1
for (file in 1:length(true_file)) {
  est_phyddle <- read.csv(paste(wkdirs, est_file[file], sep = ""))
  true_phyddle <- read.csv(paste(wkdirs, true_file[file], sep = ""))
  
  if (  grepl("fix", est_file[file], ignore.case = FALSE)) {
    fix <- 1
  } else {
    fix <-0
  }

  
  nstates <- tree_size[file] - 1 
  
  #results <- matrix(0, nrow = nstates, ncol = 3)
  #results_phyd <-matrix(0, nrow = nstates, ncol = 3)
  
  # Binary estimates
  #bin_phy <- matrix(NA, nrow = dim(true_phyddle)[1], ncol = dim(true_phyddle)[2] -1)
  
  #colnames(results_phyd) <- c("DL", "ML", "DL_ML")
  
  for (i in c(1:nstates)) {
    colZero <- 2+(i-1)*2
    colOne <- 3+(i-1)*2
    
    # There shouldn't be any NAs since the trees are the same size for a dataset 
    infer <- (est_phyddle[,colOne] == apply(est_phyddle[,colZero:colOne], 1, max)) * 1
    DL <- mean(infer == true_phyddle[,i + 1])
    
    prob_by_node[k, ] <- c(i, DL, tree_size[file], network_size[file], fix)
    k <- k + 1
  }
  
}

prob_by_node_df <- as.data.frame(prob_by_node)
prob_by_node_df$n_tips <- as.factor(prob_by_node_df$n_tips)
prob_by_node_df$network_size <- as.factor(prob_by_node_df$network_size)
prob_by_node_df$fixed<- as.factor(prob_by_node_df$fixed)

prob_by_node_df <- prob_by_node_df[- which(prob_by_node_df$n_tips == 200), ]
panelB <- ggplot(prob_by_node_df, aes(node_number, probability_correct, color = network_size, shape = fixed)) + geom_point(size = 1) + scale_shape_manual(labels = c( "variable", "fixed"), values = c(20, 5))+
  facet_wrap(~n_tips, scales = "free_x", nrow = 1) + theme_classic()+  
  labs(x = "Phyddle node number",
       y = "proportion correct", color = "network max\ntree size", shape = "network\ntree size") 
pdf("~/asr_writing/manuscript/figs/accuracy_fix_in_var.pdf", width =5, height = 3)
print(panelB)
dev.off()
png("~/asr_writing/manuscript/figs/accuracy_fix_in_var.png", width =5, height = 3, units = "in", res = 500)
print(panelB)
dev.off()

# ggplot(prob_by_node_df, aes(node_number, probability_correct, color = n_tips)) + geom_point(size = 0.5) + 
#   facet_wrap(~network_size, scales = "free_x", nrow = 1) + theme_classic()+  
#   labs(x = "Phyddle node number",
#        y = "Accuracy", color = "number of tips") 

tip_50_net_50 <- intersect(intersect (which(prob_by_node_df$n_tips == 50 ), which(prob_by_node_df$network_size == 50)), which(prob_by_node_df$fixed == 1) )
mean(prob_by_node_df$probability_correct[tip_50_net_50])
tip_50_net_50 <- intersect(intersect (which(prob_by_node_df$n_tips == 50 ), which(prob_by_node_df$network_size == 50)), which(prob_by_node_df$fixed == 0) )
mean(prob_by_node_df$probability_correct[tip_50_net_50])

tip_50_net_100 <- intersect(which(prob_by_node_df$n_tips == 50 ), which(prob_by_node_df$network_size == 100) )
mean(prob_by_node_df$probability_correct[tip_50_net_100])
tip_50_net_200 <- intersect(which(prob_by_node_df$n_tips == 50 ), which(prob_by_node_df$network_size == 200) )
mean(prob_by_node_df$probability_correct[tip_50_net_200])

tip_100_net_100 <-  intersect(intersect (which(prob_by_node_df$n_tips == 100 ), which(prob_by_node_df$network_size == 100)), which(prob_by_node_df$fixed == 1) )
mean(prob_by_node_df$probability_correct[tip_100_net_100])

tip_100_net_100 <-  intersect(intersect (which(prob_by_node_df$n_tips == 100 ), which(prob_by_node_df$network_size == 100)), which(prob_by_node_df$fixed == 0) )
mean(prob_by_node_df$probability_correct[tip_100_net_100])

tip_100_net_200 <- intersect(which(prob_by_node_df$n_tips == 100 ), which(prob_by_node_df$network_size == 200) )
mean(prob_by_node_df$probability_correct[tip_100_net_200])

tip_200_net_200 <-  intersect(intersect (which(prob_by_node_df$n_tips == 200 ), which(prob_by_node_df$network_size == 200)), which(prob_by_node_df$fixed == 1) )
mean(prob_by_node_df$probability_correct[tip_200_net_200])
tip_200_net_200 <-  intersect(intersect (which(prob_by_node_df$n_tips == 200 ), which(prob_by_node_df$network_size == 200)), which(prob_by_node_df$fixed == 0) )
mean(prob_by_node_df$probability_correct[tip_200_net_200])
