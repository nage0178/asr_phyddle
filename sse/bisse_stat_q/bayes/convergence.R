library(coda)
checkManual <- c()
#for (i in 1:10) {
for (i in 1:2500) {

  tree1 <- paste("output/", i, "_BiSSE_anc_states_results.tree", sep = "")
  if (!file.exists(tree1)) {
	  next
  }
  print(paste("run", i))
  run1 <- read.table(paste("output/", i, "_BiSSE_run_1.log", sep = ""), header = TRUE)
  colKeep <- c(10, 11, 26, 27, 28, 29)

  keep <- (floor(dim(run1)[1]*.10) + 1): dim(run1)[1]
 # mcmc1 <- as.mcmc(x = run1[keep, colKeep])
  mcmc1 <- as.mcmc(x = run1[, colKeep])
  ess1 <- effectiveSize(mcmc1)
  
  run2 <- read.table(paste("output/", i, "_BiSSE_run_2.log", sep = ""), header = TRUE)
  mcmc2 <- as.mcmc(x = run2[, colKeep])
 # mcmc2 <- as.mcmc(x = run2[keep, colKeep])

  ess2 <- effectiveSize(mcmc2)
  
  gd <- gelman.diag(list(mcmc1,mcmc2))

  if (any(ess1 < 200) || any(ess2 < 200) || any(gd$psrf[, 1] > 1.1)) {
	  checkManual <- c(checkManual, i) 
	  print(min(ess1))
	  print(min(ess2))
	  print(max(gd$psrf[,1]))
  }
}
write.csv(checkManual,"notConverged")

#print(checkManual)

