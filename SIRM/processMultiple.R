est1 <- read.csv("estimate/out.empirical_est.labels_cat.csv")
est2 <- read.csv("estimate_2/out.empirical_est.labels_cat.csv")
est3 <- read.csv("estimate_3/out.empirical_est.labels_cat.csv")
est4 <- read.csv("estimate_4/out.empirical_est.labels_cat.csv")
est5 <- read.csv("estimate_5/out.empirical_est.labels_cat.csv")

est_mean <- (est1 + est2 + est3 + est4 + est5 )/5
#print(colnames(est_mean))
#for (i in 1:99) {
#	start <- (i-1) * 5 + 2
#	end <- start + 4
#	print(sum(est_mean[1,start:end]))
#}

write.csv(est_mean, "empirical/mean_est.csv", row.names = FALSE)
