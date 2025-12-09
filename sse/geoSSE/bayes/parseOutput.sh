mkdir parseRb
for idx in {1..2500}
do

	if [ -f output/1/${idx}_ase.tre ] ; then
		echo $idx
		header=node,index
		#header=node,anc_state_1,anc_state_1_pp,anc_state_2,anc_state_2_pp 
		
		sed 's/nd/\nnd/g' ../simulate/sim.${idx}.extant.tre | grep nd | cut -f 1 -d : | cut -f 1 -d \; > nodes
	
		sed 's/)\[/\n\)\[/g' output/1/${idx}_ase.tre  | grep -F ')[' | cut -f 1 -d ","  | cut -f 2 -d = > index  #| awk -F , '{print $3 "," $4 "," $5 "," $6}' | sed -e "s/anc_state_[1..2]=//g" -e "s/anc_state_[1..2]_pp=//g" > probs
		echo $header > parseRb/node_index_${idx}
		paste -d ,  nodes index | sed 's/nd//g'  | sort -n | sed 's/^/nd/g'  >> parseRb/node_index_${idx} 
       fi	
done

# Remove- copied from elsewhere
#head -1 parseRb/node_prob_1 > all_Bayes.csv
#for idx in {1..2500}
#do
#	tail -n+2 parseRb/node_prob_${idx} >>all_Bayes.csv
#done
