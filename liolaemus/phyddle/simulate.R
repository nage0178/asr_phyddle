library(diversitree)
library(stringr)


# arguments
args        = commandArgs(trailingOnly = TRUE)
out_path    = args[1]
out_prefix  = args[2]
start_idx   = as.numeric(args[3])
batch_size  = as.numeric(args[4])
start_idx:(start_idx+batch_size-1)
rep_idx     = start_idx:(start_idx+batch_size-1)
#num_rep     = length(rep_idx)

# filesystem
tmp_fn = paste0(out_path, "/", out_prefix, ".", rep_idx)   # sim path prefix
phy_fn = paste0(tmp_fn, ".tre")               # newick file
dat_fn = paste0(tmp_fn, ".dat.csv")           # csv of data
lbl_fn = paste0(tmp_fn, ".labels.csv")        # csv of labels (e.g. params)
asr_true_fn = paste0(tmp_fn, ".anc_state.csv")   # csv of labels (e.g. params)

# Get the first state along a branch (forward time)
firstState <- function (mat) {
  mat[1,2]
}

# Get the last state along a branch (forward time)
lastState <- function(mat) {
  mat[dim(mat)[1], 2]
}

combine_states <- function(state_table) {
  paste0(state_table[1], state_table[2], state_table[3]) 
}

getHistory <- function(tree1, phy_file, sp_rm) {
  # Fix the off by one indexing for how the history is recorded 
  tree1$hist$from <- tree1$hist$from - 1
  tree1$hist$to <- tree1$hist$to - 1
  h <- history.from.sim.discrete(tree1, 0:2)
  
  # Find the first and last state for every branch
  mat_startState <- unlist(lapply(h$history, firstState))
  mat_endState   <- unlist(lapply(h$history, lastState))

  newTree <- tree1
  
  # Create a matrix to hold the names of the children at each node
  parent_child_matrix <- matrix(0, nrow = tree$Nnode, ncol = 2,
                                dimnames = list(tree$node.label, c("child1", "child2")))
  
  # For every branch (edge) in the edge matrix
  # find the children labels
  for (i in 1:nrow(tree$edge)) {
    parent_node_num <- tree$edge[i, 1]
    child_node_num <- tree$edge[i, 2]
    
    # Map node numbers to labels
    parent_label <- ifelse(parent_node_num <= length(tree$tip.label),
                           tree$tip.label[parent_node_num],
                           tree$node.label[parent_node_num - length(tree$tip.label)])
    child_label <- ifelse(child_node_num <= length(tree$tip.label),
                          tree$tip.label[child_node_num],
                          tree$node.label[child_node_num - length(tree$tip.label)])
    
    if (parent_child_matrix[parent_label, 1]  == 0) {
      parent_child_matrix[parent_label, 1] <- child_label
    } else {
      parent_child_matrix[parent_label, 2] <- child_label
      
    }
  }
  
  # Create a matrix to store the parent, children trios of states at internal nodes
  state_table <- matrix(NA, nrow = nrow(parent_child_matrix), ncol = 5)
  rownames(state_table) <- rownames(parent_child_matrix)
  colnames(state_table) <- c("parent", "child1", "child2", "child1_lb", "child2_lb")
  
  for(i in 1:nrow(parent_child_matrix)) {
    state_table[i,2] <- mat_startState[parent_child_matrix[i,1]]
    state_table[i,3] <- mat_startState[parent_child_matrix[i,2]]
    state_table[i,4] <- parent_child_matrix[i,1]
    state_table[i,5] <- parent_child_matrix[i,2]
      
    # If the node is not the root
    if (!is.na( mat_endState[rownames(parent_child_matrix)[i]])) {
      state_table[i,1] <- mat_endState[rownames(parent_child_matrix)[i]]  
      
    # If the node is the root
    } else {
      # If the two daughters are the same, the root is the
      # same state as the daughters
      if (state_table[i,2] == state_table[i,3]) {
        state_table[i,1] <- state_table[i,2]
        
        # Otherwise, the root is widespread
        # This is only for 2 regions
      } else {
        state_table[i,1] <- 0
      }
    }
    
  }

  num_tips <- Ntip(tree)
  num_nodes <- tree$Nnode
  
  # Select half of the nodes in the tree to rotate
  # This is needed for phyddle to have examples of both orderings
  nodes_to_rotate <- sample((num_tips + 1):(num_tips + num_nodes),
                            size = floor(num_nodes * 0.5),
                            replace = FALSE)
  
  for (node in nodes_to_rotate) {
    
    # Figure out which branches to rotate based on the parent index
    branchIndex <- which(newTree$edge[,1] == node)
    one <- branchIndex[1]
    two <- branchIndex[2]
    
    # Rotate branches
    tmp <-newTree$edge[one, ]
    newTree$edge[one, ] <- newTree$edge[two, ]
    newTree$edge[two, ] <- tmp
    
    # Rotate lengths
    tmp <- newTree$edge.length[one]
    newTree$edge.length[one] <- newTree$edge.length[two]
    newTree$edge.length[two] <- tmp
    
    # Update the table of left/right descendants 
    nodeName <- newTree$node.label[node-num_tips]
    tmp <- state_table[nodeName , 2]
    state_table[nodeName, 2] <- state_table[nodeName , 3]
    state_table[nodeName, 3] <- tmp
    
    tmp <- state_table[nodeName , 4]
    state_table[nodeName, 4] <- state_table[nodeName , 5]
    state_table[nodeName, 5] <- tmp
  }
  
  comb_state <- apply(state_table, 1, combine_states)
  newStates <- sapply(comb_state, switch, 
        "111" = "1", 
        "222" = "2", 
        "012" = "3",
        "021" = "4",
        "001" = "5",
        "010" = "6",
        "002" = "7",
        "020" = "8")

  # The way phyddle is working right now the order of these things matter
  # In the current hard coding, 1 is range A , 2 is range B, 3 is widespread
  # 's/,1,1,1$/,1/g'  
  # 's/,2,2,2$/,2/g' 
  # 's/,3,1,2$/,3/g' 
  # 's/,3,2,1$/,4/g'  
  # 's/,3,3,1$/,5/g' 
  # 's/,3,1,3$/,6/g' 
  # 's/,3,3,2$/,7/g' 
  # 's/,3,2,3$/,8/g'
  
  #names <- paste0("sp", sp_rm)
  sub_tree <- newTree
  #write.tree(sub_tree, file = "tmp.tre")

  #sub_tree <- drop.tip(newTree, sp_rm)
  for(i in sp_rm) {
    sub_tree <- drop.tip(sub_tree, i)
  }

  write.tree(sub_tree, file = phy_file)

  return (cbind(newStates, state_table[,4:5]))
}


# Labels from Goldberg GeoSSE paper
label_names = c("sC", "sF", "sCF", "xC", "xF", "dC", "dF", "frac_C", "frac_F", "frac_CF")
#pars is ordered sA, sB, sAB, xA, xB, dA, dB.

nTrees <- start_idx - 1
nullTrees <- 0
noVarTree <- 0
set.seed(nTrees)
totTrees <- start_idx+batch_size-1
params <- matrix (0, nrow = totTrees, ncol = 7)
i <- 1
while (nTrees < totTrees) {
  
  # Priors  are from Goldberg GeoSSE paper
  sC <- rexp(1, 1/.1)
  sF <- rexp(1, 1/.1)
  xC <- rexp(1, 1/.1)
  xF <- rexp(1, 1/.1)
  dC <- rexp(1, 1/.1)
  dF <- rexp(1, 1/.1)
  sCF <- rexp(1, 1/.1)
  
  pars <- c(sC, sF, sCF, xC, xF, dC, dF)
  
  # For Cean clade 
  tree_tips <- 52
  max_tax <- tree_tips
  # sampled proportions of C, F, and CF: 0.89, 0.95, and 0.91
#  samp_frac  <- runif(1, .85, 1)
#  max_tax <- round(tree_tips / samp_frac)
  
  
  tree <- tree.geosse(pars, max.taxa=max_tax, max.t=Inf, include.extinct=TRUE,
              x0=NA)
  extant <-  which(str_detect(tree$tip.label, "sp"))

  # Check if tree exists
  if (is.null(tree) | sum(str_detect(tree$tip.label, "sp"))  != max_tax ) {
    nullTrees <- nullTrees + 1
    
  # Confirm at least two tip states present
  } else if (length(table(tree$tip.state[extant])) < 3) {
    noVarTree <- noVarTree + 1
    
    # Write all relevant files for phyddle
    } else {

    # Select extant tips to remove
    sp_rm <- sample(tree$tip.label[which(str_detect(tree$tip.label, "sp"))], 
         size = max_tax - tree_tips,
         replace = FALSE)

    stateHistory <- getHistory(tree, phy_fn[i], sp_rm)
    # This file should actually be in the form
    # node, state(1-8), node l, node r
    write.table(stateHistory, file = asr_true_fn[i],quote = FALSE, col.names = FALSE, sep = "," )  
    
    
    # Change this
    df_state = data.frame(taxa=tree$tip.label, data=tree$tip.state, region1 = 0, region2 = 0)
 #   df_state = data.frame(taxa=tree$tip.label, data = sapply(tree$tip.state, switch, 
 #  			"0" = "1,1",
 #  			"1" = "1,0",
 # 			"2" = "0,1"))
    # Calculate sampling fractions by type before you remove
    #df_state <- df_state[-row_rm, ]
    
    df_state$region1[df_state$data == 0 | df_state$data == 1] <- 1
    df_state$region2[df_state$data == 0 | df_state$data == 2] <- 1
    
    extant <- which(str_detect(df_state$taxa, "sp"))
    df_state <- df_state[extant, ]

    row_rm <- which(df_state$taxa %in% sp_rm)
    #widespread
    frac_0 <- 1
    # in 1 
    frac_1 <- 1
    # in 2
    frac_2 <- 1
    
    # These are going to include extinct species- need to recalculate
    if (length(row_rm) > 0) {
      remain <- df_state[-row_rm, ]
      frac_0 <- sum(remain$region1 + remain$region2 == 2)/ sum(df_state$region1 + df_state$region2 == 2) 
      frac_1 <- sum(remain$region1 == 1 & remain$region2 == 0)/ sum(df_state$region1 == 1 &  df_state$region2 == 0) 
      frac_2 <- sum(remain$region1 == 0 & remain$region2 == 1)/ sum(df_state$region1 == 0 &  df_state$region2 == 1) 
      df_state <- df_state[-row_rm, ]
      #print(sum(df_state$region1 + df_state$region2 == 2))
      #print(sum(df_state$region1 == 1 &  df_state$region2 == 0))
      #print(sum(df_state$region1 == 0 &  df_state$region2 == 1))
    } 
    
    df_state <- df_state[ ,-2]

    write.csv(df_state, file=dat_fn[i], row.names=F, quote=F)
      
    # save learned labels (e.g. estimated data-generating parameters)
    pars <- c(log(pars, base =10), frac_1, frac_2, frac_0)
    names(pars) = label_names
    df_label = data.frame(t(pars))
    write.csv(df_label, file=lbl_fn[i], row.names=F, quote=F)
    
    # This is to look at distributions of parameters after
    # rejection in R
    #params[nTrees + 1, ] <- pars
    nTrees <- nTrees + 1
    
    # Want to same trees independent of how phyddle is run
    set.seed(nTrees)
    i <- i+1
    
  }
}

