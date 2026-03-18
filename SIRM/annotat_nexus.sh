idx=1
nstate=5

outfile=empirical/out.${idx}.ann.nex
infile=empirical/out.${idx}.form.tre
nodeLabel=empirical/out.${idx}.node_labels.csv
dat=empirical/out.${idx}.dat.csv
est=empirical/mean_est.csv

# Beginning of nexus format
echo \#NEXUS > $outfile 
echo Begin taxa\; >> $outfile 
taxa=$(sed 's/(/\n/g' $infile | wc | awk '{print $1}')
echo Dimensions ntax = ${taxa}\; >> $outfile 

# Write taxa names for extant only tree
echo Taxlabels >> $outfile
sed -e 's/:/:\n/g' -e 's/,/,\n/g' -e 's/)/)\n/g' -e 's/(/(\n/g' $infile | grep :  | grep -v Node | awk -F : '{print $1}' >> $outfile

echo \; >> $outfile 
echo End\; >> $outfile 


# Begin tree section
echo Begin trees\; >> $outfile

# Start with the extant only tree
newLine=$(cat $infile)  
echo tree TREE1 = $newLine  >> $outfile

((index=taxa+1))

# Create annotation for each internal node
# This assumes that all the nodes have probability for all possible states
for line in $(tail -n+2 $nodeLabel) 
do

	# Names in the tree will be replaced with names from phyddle
	nodeOld=$(echo $line | awk -F , '{print $1}')
	nodeNew=$(echo $line | awk -F , '{print $2}')
	
	stateString="index=${index}"
	probString=
	
	i=1
	((col=(nodeNew)*nstate+1+i))

	# Order the probabilities, print the state number next to the probabiliites 
	ordered=$(grep ^${idx}, $est | awk -F , -v var1=$col -v var2=$nstate '{ for (i = 0; i < var2; i++) printf("%s,%d\n", $(i+var1), i)}' | sort -g -r)
#	echo $ordered

	# Combines annoations for one node
	for line in $(echo $ordered |cut -d ' ' -f 1-3) 
	do
		
		prob=$(echo $line | awk -F ,  '{print $1}')
		state=$(echo $line | awk -F ,  '{print $2}')
		stateString="${stateString},anc_state_${i}=${state}"
		probString="${probString},anc_state_${i}_pp=${prob}"
		((i=i+1))
		
	done

	other=$(echo $ordered | sed 's/ /\n/g' | tail -n+4 | awk -F',' '{sum+=$1} END{print sum;}')

	# This will only work for up to three states
	probString="${probString},anc_state_other_pp=${other}"

	sed -i -e "s/)${nodeOld}:/)[${stateString}${probString}]:/" -e "s/)${nodeOld};/)[${stateString}${probString}];/" $outfile
	
	((index=index+1))
done

# Create tip annotations
# This will depend on the format of the input data
# The problem is your file is not in csv format
index=1
#for line in $(cat ${dat}) 
#do
while IFS= read -r line; do

	tip=$(echo $line| awk '{print $1}')
	state=$(echo $line| awk '{print $2}')

	state_1=NA
	state_2=NA
	state_3=NA
	prob1=1.0
	prob2=0.0
	prob3=0.0
	
	if [[ "$state" ==  "10000"  ]]; then
		state_1=0
	fi

	if [[ "$state" ==  "01000"  ]]; then
		state_1=1
	fi

	if [[ "$state" ==  "00100"  ]]; then
		state_1=2
	fi

	if [[ "$state" ==  "00010"  ]]; then
		state_1=3
	fi

	if [[ "$state" ==  "00001"  ]]; then
		state_1=4
	fi

	echo $state_1
	
	
	var="[index=${index},anc_state_1=${state_1},anc_state_2=${state_2},anc_state_3=${state_3},anc_state_1_pp=${prob1},anc_state_2_pp=${prob2},anc_state_3_pp=${prob3},anc_state_other_pp=0.0]"
	((index=index+1))

	# Combine these
	sed -i -e "s/,${tip}:/,${tip}${var}:/" -e "s/(${tip}:/(${tip}${var}:/" $outfile


done < $dat
echo End\; >> $outfile 
sed -i 's/index/\&index/g' $outfile
