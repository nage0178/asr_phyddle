library(coda)
dirs <- c("fix/50/bayes/", 
          "fix/100/bayes/",
          "fix/200/bayes/",
          "var/50/bayes/",
          "var/100/bayes/",
          "var/200/bayes/"
	  )


for (dir in dirs) {
	print(dir)
	checkManual <- c()
	#for (i in 1:10) {
	for (i in 1:2500) {
	  run1 <- read.table(paste(dir, "output/", i, "_run_1.log", sep = ""), header = TRUE)
	  colRm <- c(1, 4, 6, 7)

  	  keep <- (floor(dim(run1)[1]*.25) + 1): dim(run1)[1]
	  mcmc1 <- as.mcmc(x = run1[keep, -colRm])
	  ess1 <- effectiveSize(mcmc1)
	  
	  run2 <- read.table(paste(dir, "output/", i, "_run_2.log", sep = ""), header = TRUE)
	  mcmc2 <- as.mcmc(x = run2[keep, -colRm])
	
	  ess2 <- effectiveSize(mcmc2)
	  
	  gd <- gelman.diag(list(mcmc1,mcmc2))
	
	  if (any(ess1 < 200) || any(ess2 < 200) || any(gd$psrf[, 1] > 1.1)) {
		  checkManual <- c(checkManual, i) 
	  }
	}
	write.csv(checkManual, paste(dir, "notConverged", sep = ""))

}
#print(checkManual)

