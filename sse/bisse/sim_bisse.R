#!/usr/bin/env Rscript
#library(phytools)
library(castor)
library(ape)
library(diversitree)
library(R.utils)

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

# dataset setup
num_states = 2
label_names = c( paste0("log_birth_",1:num_states), "log_death_1", "log_death_2", "log_state_rate_1", "log_state_rate_2",  "start_state")

num_taxa <-  0
max_taxa <-  50
numStarts <- 25

jitter.pars <- function(pars) {
	multfac <- rep(0, 6) 

        log_birth = runif(n=num_states, -2, 0)
        multfac[1] = 10^log_birth

        log_birth = runif(n=num_states, -2, 0)
        multfac[2] = 10^log_birth

        multfac[3]= min(multfac[1]) * 10^runif(n=1, -2, 0)
        multfac[4]= min(multfac[2]) * 10^runif(n=1, -2, 0)

        log_state_rate = runif(1,-2, 0)
        multfac[5]= 10^log_state_rate

        log_state_rate = runif(1,-2, 0)
        multfac[6]= 10^log_state_rate

	return(multfac)
}

# simulate each replicate
for (i in 1:num_rep) {

	# set RNG seed
	set.seed(rep_idx[i])
	print(rep_idx[i])

	# Do not resimulate the file exists
	if (file.exists(asr_fn[i]) ) {
	        print(asr_fn[i])
	        next
	}
	
	num_taxa <- -1
	# The tree is not the correct size
	while (num_taxa !=  max_taxa) {
	
	    # simulate parameters
	    start_state = sample(1:2, size=1)
	    log_state_rate_1 = runif(1,-3, 0)
	    log_state_rate_2 = runif(1,-3, 0)
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
	    print(Q)
	
	    log_birth = runif(n=num_states, -2, 0)
	    birth = 10^log_birth
	
	    death = min(birth) * 10^runif(n=1, -2, 0)
	    death = rep(death, num_states)
	    log_death = log(death[1], base=10)
	    parameters = list(
	        birth_rates=birth,
	        death_rates=death,
	        transition_matrix_A=Q
	    )
	
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
	
	# Infer ancestral states
	if (rep_idx[i] <= 2500) {
	    print(rep_idx[i])
	    print("calculating MLE")
    
    		if (mean(state_sim) == 0 ||  mean(state_sim) == 1) {
    		        asr_est_print <- matrix(state_sim[-1])

    		} else {
    			print(pars)
    			lik.s <- make.bisse(res_sim$tree, state_sim)
    			fit.s <- find.mle(lik.s, pars, method = "optim",  lower = 0, upper = 10, control = list( reltol = 1e-6))
    			best.fit.s <- fit.s 
    			best.lnL <- best.fit.s$lnLik
    			fitParams <- array(c(best.lnL, best.fit.s$par))

    		    	newPars <-pars
    			for (j in 1:numStarts) {
    			        print("Iteration")
    			        print(j)
    		    		print(format(fit.s$par, scientific = F ))
    		    		newPars <- jitter.pars(pars)

				fit.s <- find.mle(lik.s, newPars, method = "optim", lower = 0, upper = 10, control = list( reltol = 1e-6))
    			        fitParms <- rbind(fitParams, c(lnL =fit.s$lnLik, fit.s$par))

    			        if (fit.s$lnLik > best.lnL && !any(fit.s$par > 10)) {
    			    	    print("improvement")
    			    	    best.fit.s <- fit.s
    			    	    best.lnL <- fit.s$lnLik

    			        }
    			}

    			st.s <- t(asr.marginal(lik.s, coef(best.fit.s)))
    			asr_est_print <- data.frame(cbind((st.s[, 2] > .5) * 1, st.s))
    			#asr_est_print <- data.frame((asr_est$ancestral_states - 1))
    			#
    			print("final")
    			print(best.fit.s$par)
    		}
    		rownames(asr_est_print) <- tree_sim$node.label
    		write.table(asr_est_print, file = asr_fn[i], sep = ',', col.names =F, row.names =T, quote = F)
    	}
}

# done!
quit()
