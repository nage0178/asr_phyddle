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
	print(mean(root_est[1:2500] == truth[1:2500, 2]))
}

