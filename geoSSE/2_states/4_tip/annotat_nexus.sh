idx=$1
nstate=3

# Beginning of nexus format
echo \#NEXUS >  simulate/sim.${idx}.ann.nex
echo Begin taxa\; >> simulate/sim.${idx}.ann.nex
taxa=$(sed 's/(/\n/g' simulate/sim.${idx}.extant.tre | wc | awk '{print $1}')
echo Dimensions ntax = ${taxa}\; >> simulate/sim.${idx}.ann.nex

# Write taxa names for extant only tree
echo Taxlabels >> simulate/sim.${idx}.ann.nex
sed -e 's/:/:\n/g' -e 's/,/,\n/g' -e 's/)/)\n/g' -e 's/(/(\n/g' simulate/sim.1.extant.tre | grep :  | grep -v node | awk -F : '{print $1}' >> simulate/sim.${idx}.ann.nex

echo \; >> simulate/sim.${idx}.ann.nex
echo End\; >> simulate/sim.${idx}.ann.nex


# Begin tree section
echo Begin trees\; >> simulate/sim.${idx}.ann.nex

# Start with the extant only tree
newLine=$(cat simulate/sim.${idx}.extant.tre)  
echo tree TREE1 = $newLine  >> simulate/sim.${idx}.ann.nex

((index=taxa+1))

# Create annotation for each internal node
for line in $(tail -n+2  simulate/sim.${idx}.node_labels.csv) 
do

	# Names in the tree will be replaced with names from phyddle
	nodeOld=$(echo $line | awk -F , '{print $1}')
	nodeNew=$(echo $line | awk -F , '{print $2}')
	
	stateString="index=${index}"
	probString=
	
	i=1
	((col=(nodeNew)*nstate+1+i))

	# Order the probabilities, print the state number next to the probabiliites 
	ordered=$(grep ^${idx}, estimate/sim.test_est.labels_cat.csv | awk -F , -v var1=$col -v var2=$nstate '{ for (i = 0; i < var2; i++) printf("%s,%d\n", $(i+var1), i)}' | sort -n -r)
	
	# Combines annoations for one node
	for line in $(echo $ordered) 
	do
		
		prob=$(echo $line | awk -F ,  '{print $1}')
		state=$(echo $line | awk -F ,  '{print $2}')
		stateString="${stateString},anc_state_${i}=${state}"
		probString="${probString},anc_state_${i}_pp=${prob}"
		((i=i+1))
		
	done
	# This will only work for up to three states
	probString="${probString},anc_state_other_pp=0.0"

	sed -i "s/)${nodeOld}:/)[${stateString}${probString}]:/" simulate/sim.${idx}.ann.nex
	
	
	((index=index+1))
	echo
done

#sed -i 's/,]/]/g' simulate/sim.${idx}.ann.nex

# Create tip annotations
# This will depend on the format of the input data
index=1
for line in $(cat simulate/sim.${idx}.dat.csv) 
do
	tip=$(echo $line| awk -F , '{print $1}')
	region_1=$(echo $line| awk -F , '{print $2}')
	region_2=$(echo $line| awk -F , '{print $3}')
	state_1=NA
	state_2=NA
	state_3=NA
	prob1=1.0
	prob2=0.0
	prob3=0.0
	
	if [[ "$region_1" ==  1  && "$region_2" == 0 ]]; then
		state_1=0
	fi

	if [[ "$region_1" ==  1  && "$region_2" == 1 ]]; then
		state_1=2
	fi

	if [[ "$region_1" ==  0  && "$region_2" == 1 ]]; then
		state_1=1
	fi
	
	var="[index=${index},anc_state_1=${state_1},anc_state_2=${state_2},anc_state_3=${state_3},anc_state_1_pp=${prob1},anc_state_2_pp=${prob2},anc_state_3_pp=${prob3},anc_state_other_pp=0.0]"
	((index=index+1))

	# Combine these
	sed -i -e "s/,${tip}:/,${tip}${var}:/" -e "s/(${tip}:/(${tip}${var}:/" simulate/sim.${idx}.ann.nex


done
echo End\; >> simulate/sim.${idx}.ann.nex
sed -i 's/index/\&index/g' simulate/sim.${idx}.ann.nex
