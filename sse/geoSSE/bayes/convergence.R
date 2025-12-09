library(coda)
checkManual <- c()
#for (i in 1:10) {
for (i in 1:2500) {
  file1 <- paste("output/1/", i, "_model.log", sep = "")
  file2 <- paste("output/2/", i, "_model.log", sep = "")
  
  tree1 <- paste("output/1/", i, "_ase.tre", sep = "")
  tree2 <- paste("output/2/", i, "_ase.tre", sep = "")

  if (file.exists(tree1) && file.exists(tree2)) {
  	print(i)

  	run1 <- read.table(file1, header = TRUE)
  	colKeep <- c(5, 6, 13, 14, 15, 26, 27)

  	keep <- (floor(dim(run1)[1]*.25) + 1): dim(run1)[1]
  	mcmc1 <- as.mcmc(x = run1[keep, colKeep])
  	ess1 <- effectiveSize(mcmc1)
  	
  	run2 <- read.table(file2, header = TRUE)
  	mcmc2 <- as.mcmc(x = run2[keep, colKeep])

  	ess2 <- effectiveSize(mcmc2)
  	
  	gd <- gelman.diag(list(mcmc1,mcmc2))

  	if (any(ess1 < 200) || any(ess2 < 200) || any(gd$psrf[, 1] > 1.1)) {
  	        checkManual <- c(checkManual, i) 
        	print(min(ess1))
        	print(min(ess2))
        	print(max(gd$psrf[, 1]))
  	}
  }
}
write.csv(checkManual,"notConverged")

print(checkManual)

