for dir in {fix,var}/{50,100,200}/bayes
do

	cd $dir
	size=$(echo $dir | cut -d '/' -f 2 ) 

	mkdir parseRb

	for idx in {1..2500}
	do
		echo $idx
		header=node,anc_state_1,anc_state_1_pp,anc_state_2,anc_state_2_pp 
		
		# Make list of all internal nodes in the order they appear in the tree
		sed 's/node/\nnode/g' ../simulate/out.${idx}.tre | grep node | cut -f 1 -d : | cut -f 1 -d \; > nodes

		# Make list of all node probabilities in the order they appear in the tree
		sed 's/)\[/\n\)\[/g' output/${idx}_ase.tree  | grep -F ')[' | awk -F , '{print $3 "," $4 "," $5 "," $6}' | sed -e "s/anc_state_[1..2]=//g" -e "s/anc_state_[1..2]_pp=//g" > probs
		echo $header > parseRb/node_prob_${idx}

		# Sort the nodes based on their index
		paste -d ,  nodes probs  | sed 's/node//g' | sort -n | sed 's/^/node/g'  >> parseRb/node_prob_${idx} 

		while [[ $(wc -l parseRb/node_prob_${idx} | cut -d " " -f 1 ) -lt $size ]]; do
		  echo "node,,,," >>  parseRb/node_prob_${idx}
		done
	done
	
	head -1 parseRb/node_prob_1 > all_Bayes.csv
	for idx in {1..2500}
	do
		tail -n+2 parseRb/node_prob_${idx} >> all_Bayes.csv
	done
	rm nodes
	rm probs

	cd ../../../
	#exit

done
