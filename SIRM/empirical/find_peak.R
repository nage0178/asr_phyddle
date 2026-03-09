library(phytools)

dat <- as.data.frame(read.csv("EBOV_maxCases.csv"))
dat_SLE <- dat[which(dat$country == "SLE"), ]

# Add the districts together that correspond to the same regions in the analysis
pop1_row <- c(1, 4, 6)
pop2_row <- c(8, 7, 2)
pop3_row <- c(12, 10, 2)
pop4_row <- c(11, 9, 3)
pop5_row <- c(13, 14)

comb_pop <- matrix(0, nrow = 5, ncol = 99)

for (i in 1:length(pop1_row)) {
  comb_pop[1, ] <- comb_pop[1, ] + unlist(dat_SLE[pop1_row[i], 4:102])
  comb_pop[2, ] <- comb_pop[2, ] + unlist(dat_SLE[pop2_row[i], 4:102])
  comb_pop[3, ] <- comb_pop[3, ] + unlist(dat_SLE[pop3_row[i], 4:102])
  comb_pop[4, ] <- comb_pop[4, ] + unlist(dat_SLE[pop4_row[i], 4:102])

    if (i < 3) {
      comb_pop[5, ] <- comb_pop[5, ] + unlist(dat_SLE[pop5_row[i], 4:102])
    }
}

# First cases are on week 20 
peak_incidence <- apply(comb_pop, 1, function(x){which(x == max(x))}) - 20 


prevalence_1week <- matrix(0, nrow = dim(comb_pop)[1], ncol = dim(comb_pop)[2])
prevalence_2week <- matrix(0, nrow = dim(comb_pop)[1], ncol = dim(comb_pop)[2])

for (i in 3:dim(comb_pop)[2]) {
  # Note the first two weeks are zeros anyway
  prevalence_1week[, i] <- comb_pop[, i] + comb_pop[, i - 1]
  prevalence_2week[, i] <- comb_pop[, i] + comb_pop[, i - 1] + comb_pop[, i -2]
  
}

# First cases are on week 20 
peak_prev_1week <- apply(prevalence_1week, 1, function(x){which(x == max(x))}) - 20 
peak_prev_2week <- apply(prevalence_2week, 1, function(x){which(x == max(x))}) - 20 

