

# This gets the parent daughter relationships in the tree 
# using the indeces from rb

mkdir -p parseOutput
# Rm the rb part
RB_CMD=~/revbayes/projects/cmake/rb

#for idx in {1..2500}
#do
#	RB_STR="IDX=${idx}; source(\"scripts/parse.Rev\")"
#	echo ${RB_STR} | ${RB_CMD} 
#done

# This creates a file with that matches the indeces from rb
# and the internal node names from the original tree
./parseOutput.sh

# Gets the point estimate (out of the 8 possible parent/daughter
# triplets) based on the most common frequent from the MCMC for 
# each replicate
Rscript process_anc_state.R

# Makes a confusion matrix with the point estimates
#Rscript makeConfusion.R 
#
Rscript get_probs.R
#
Rscript compareProbs.R
