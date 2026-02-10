#!/usr/bin/env Rscript
#library(phytools)
library(castor)
library(ape)
library(diversitree)
library(R.utils)
library(dispRity)

# disable warnings
#options(warn = -1)

# arguments
args        = commandArgs(trailingOnly = TRUE)
out_path    = args[1]
out_prefix  = args[2]
start_idx   = as.numeric(args[3])
batch_size  = as.numeric(args[4])
rep_idx     = start_idx:(start_idx+batch_size-1)
num_rep     = length(rep_idx)
get_mle     = TRUE

# filesystem
tmp_fn = paste0(out_path, "/", out_prefix, ".", rep_idx)   # sim path prefix
phy_fn = paste0(tmp_fn, ".tre")               # newick file
dat_fn = paste0(tmp_fn, ".dat.csv")           # csv of data
lbl_fn = paste0(tmp_fn, ".labels.csv")        # csv of labels (e.g. params)
asr_true_fn = paste0(tmp_fn, ".anc_state.csv")     # csv of labels (e.g. params)
asr_fn = paste0(tmp_fn, ".asr.csv")        # csv of labels (e.g. params)
age_prop_fn = paste0(tmp_fn, ".age_prop.csv")        # csv of node ages 

# dataset setup
num_states = 2
label_names = c( paste0("log_birth_",1:num_states), "log_death_1", "log_death_2", "log_state_rate_1", "log_state_rate_2",  "start_state")

num_taxa <-  0
max_taxa <-  50
numStarts <- 25

find_stationary <- function (birth, death, state) {

 g <- birth[1] - death[1] - birth[2] + death[2]
 q01 <- state[1]
 q10 <- state[2]
 a <- -g 
 b <- g - q01 -q10
 c <- q10
 root1 <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
 root2 <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
 root <- c(root1, root2)


 root[intersect(which(root <=1) , which(root >=0))]

}


# simulate each replicate
for (i in 1:num_rep) {

	# set RNG seed
	set.seed(rep_idx[i])
	print(rep_idx[i])

	
	num_taxa <- -1
	# The tree is not the correct size
	while (num_taxa !=  max_taxa) {
	
	    # simulate parameters
	    log_state_rate_1 = runif(1,-2, 0)
	    log_state_rate_2 = runif(1,-2, 0)
	    #log_state_rate_1 = runif(1,-3,-0.2)
	    #log_state_rate_2 = runif(1,-3,-0.2)
	    #log_state_rate = runif(1,-3,0)

	    state_rate_1 = 10^log_state_rate_1
	    state_rate_2 = 10^log_state_rate_2

	    Q = matrix(state_rate_1,
	               ncol=num_states, nrow=num_states)
	    Q[2,1] <- state_rate_2
	    diag(Q) = 0
	    diag(Q) = -rowSums(Q)
	
	    log_birth = runif(n=num_states, -2, 0)
	    birth = 10^log_birth
	
	    death = birth * 10^runif(n=2, -2, 0)
	    #death = rep(death, num_states)
	    log_death = log(death[1], base=10)
	    parameters = list(
	        birth_rates=birth,
	        death_rates=death,
	        transition_matrix_A=Q
	    )
	    station_0 <- find_stationary(birth, death, c(state_rate_1, state_rate_2))
	    draw <-  runif(1)
	    start_state <- (draw > station_0) * 1 + 1
	
	    # simulate tree/data
	    res_sim = simulate_dsse(
	            Nstates=num_states,
	            parameters=parameters,
	            start_state=start_state,
	            sampling_fractions=1,
	            max_extant_tips=max_taxa,
	            include_labels=T,
	            no_full_extinction=F)
	
	    # check if tree is valid
	    num_taxa = length(res_sim$tree$tip.label)
	}
	
	# save tree
	tree_sim <- makeNodeLabel(res_sim$tree, method = "number", prefix = "node")
	write.tree(tree_sim, file=phy_fn[i])
	
	# save data
	# save learned labels (e.g. estimated data-generating parameters)
	
	# This is saving the tip states
	state_sim = res_sim$tip_states - 1
	state_sim_node = res_sim$node_states - 1
	df_state = data.frame(taxa=tree_sim$tip.label, data=state_sim)
	write.csv(df_state, file=dat_fn[i], row.names=F, quote=F)
	
	# These are the simulation parameters
	label_sim = c( log_birth[1], log_birth[2], log_death[1], log_death[2], log_state_rate_1, log_state_rate_2,  start_state-1)
	
	names(label_sim) = label_names
	df_label = data.frame(t(label_sim))
	write.csv(df_label, file=lbl_fn[i], row.names=F, quote=F)
	
	# These are the ancestral states
	df_label = data.frame(t(label_sim))
	anc_state <- data.frame(t(res_sim$node_states -1 ))
	
	colnames(anc_state) <- paste("node", c(1:(max_taxa -1)), sep = "")
	write.table(t(anc_state), file = asr_true_fn[i], quote =F, col.names =F, sep = ',')
	
	pars <- c(birth[1], 
		  birth[2],
	    	  death[1],  
	    	  death[2],  
	    	  state_rate_1,  
	    	  state_rate_2
	    	 )

	ages <- (tree.age(tree_sim))
	ages <- ages[(max_taxa+1):dim(ages)[1],]
	#ages <- ages[(dim(asr_est_print)[1]+2):dim(ages)[1],]

	agesProp <- ages[,1] / max(ages[,1])
	names(agesProp) <- ages[,2]
	agesProp <- t(agesProp)
	write.csv(agesProp, file = age_prop_fn[i], quote= F, row.names =F)
	
}

warnings()
# done!
quit()
