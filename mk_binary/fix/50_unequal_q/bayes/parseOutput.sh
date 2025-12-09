mkdir parseRb
for idx in {1..2500}
do
	echo $idx
	header=node,anc_state_1,anc_state_1_pp,anc_state_2,anc_state_2_pp 
	
	sed 's/node/\nnode/g' ../simulate/out.${idx}.tre | grep node | cut -f 1 -d : | cut -f 1 -d \; > nodes
	sed 's/)\[/\n\)\[/g' output/${idx}_ase.tree  | grep -F ')[' | awk -F , '{print $3 "," $4 "," $5 "," $6}' | sed -e "s/anc_state_[1..2]=//g" -e "s/anc_state_[1..2]_pp=//g" > probs
	echo $header > parseRb/node_prob_${idx}
	paste -d ,  nodes probs  | sed 's/node//g' | sort -n | sed 's/^/node/g'  >> parseRb/node_prob_${idx} 
done

head -1 parseRb/node_prob_1 > all_Bayes.csv
for idx in {1..2500}
do
	tail -n+2 parseRb/node_prob_${idx} >>all_Bayes.csv
done
