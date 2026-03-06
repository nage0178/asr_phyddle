setwd("~/asr_phyddle/mk_binary/fix/")
for (i in c(50,100,200)) {
	print(i)
	print("root alone")
	est   <- read.csv(paste("root_alone/", i, "/estimate/out.test_est.labels_cat.csv" , sep = ""))
	truth <- read.csv(paste("root_alone/", i, "/estimate/out.test_true.labels_cat.csv", sep = ""))
	root_est <-  1 * (est[, 3] > .5)
	print(mean(root_est == truth[, 2]))
	
	print("all nodes")
	est   <- read.csv(paste(i, "/est.csv"  , sep = ""))
	truth <- read.csv(paste(i, "/truth.csv", sep = ""))
	root_est <-  1 * (est[, 3] > .5)
	print(mean(root_est == truth[, 2]))
}

#est <- read.csv("~/asr_phyddle/mk_binary/fix/100_root/est.csv")
#truth <- read.csv("~/asr_phyddle/mk_binary/fix/100_root/truth.csv")
#
## This should be the root based on the node number, this is the re-ordered dataset
#root_est <-  1 * (est[, 3] > .5)
#mean(root_est == truth[, 2])
#
phy_true <- read.csv("~/asr_phyddle/mk_binary/fix/100_root/estimate/out.test_true.labels_cat.csv")
phy_est <- read.csv("~/asr_phyddle/mk_binary/fix/100_root/estimate/out.test_est.labels_cat.csv")
# Col 1 is idx, col 2 is state 0 for root, col3, is state 1 for root
phy_root_est <- 1 * (phy_est[, 3] > .5 )
mean(phy_root_est == phy_true[, 2])
#
#phy_true <- read.csv("~/asr_phyddle/mk_binary/fix/50_root/estimate/out.test_true.labels_cat.csv")
#phy_est <- read.csv("~/asr_phyddle/mk_binary/fix/50_root/estimate/out.test_est.labels_cat.csv")
## Col 1 is idx, col 2 is state 0 for root, col3, is state 1 for root
#phy_root_est <- 1 * (phy_est[, 3] > .5 )
#mean((phy_root_est == phy_true[, 2])[1:2500])
