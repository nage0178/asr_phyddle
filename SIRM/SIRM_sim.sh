python3 post_exponentialphase_SIRM_sim_one.py $1 $2 $3 $4 #>> tmp 
wait

((total=$3+$4))
for ((i=$3;i<$total;i++))
do

	# Move the tree so you don't overwrite it 
	mv $1/$2.$i.tre $1/$2.$i.org.tre

	# Reformat labels
        #./labels_one.sh $i & 

	# Label the internal nodes
	python3 label.py $1/$2.$i.nex.tre $1/$2.$i.labeled.tre $1/$2.$i.tmp_anc_state.csv $1/$2.$i.tre 

	# Remove the extra parenthesis for the tips
	#./rmExtraParen.sh $1/$2.$i.labeled.tre $1/$2.$i.nex.tre
	sed 's/Node/\nNode/g' $1/$2.$i.tre | sed "s/':/,\n/g" | grep Node > $1/nodes_$i
	grep -f $1/nodes_$i $1/$2.$i.tmp_anc_state.csv > $1/$2.$i.anc_state.csv
	rm $1/nodes_$i 


done
