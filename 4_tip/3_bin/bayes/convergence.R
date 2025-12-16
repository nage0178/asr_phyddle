library(coda)
checkManual <- c()
for (i in 1:2500) {
  run1 <- read.table(paste("output/", i, "_run_1.log", sep = ""), header = TRUE)
  colRm <- c(1, 4, 6, 7)

  mcmc1 <- as.mcmc(x = run1[, -colRm])
  ess1 <- effectiveSize(mcmc1)
  
  run2 <- read.table(paste("output/", i, "_run_2.log", sep = ""), header = TRUE)
  mcmc2 <- as.mcmc(x = run2[, -colRm])

  ess2 <- effectiveSize(mcmc2)
  
  gd <- gelman.diag(list(mcmc1,mcmc2))

  if (any(ess1 < 200) || any(ess2 < 200) || any(gd$psrf[, 1] > 1.1)) {
	  print(i)
	  checkManual <- c(checkManual, i) 
  }
}

warnings()

print(checkManual)
