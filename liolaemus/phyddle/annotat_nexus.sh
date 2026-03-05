idx=1
nstate=8

outfile=empirical/lio.${idx}.ann.nex
infile=empirical/lio.${idx}.form.tre
nodeLabel=empirical/lio.${idx}.node_labels.csv
dat=empirical/lio.${idx}.dat.csv
est=estimate/sim.empirical_est.labels_cat.csv

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
	echo $ordered

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
index=1
for line in $(cat ${dat}) 
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
	sed -i -e "s/,'${tip}':/,'${tip}'${var}:/" -e "s/('${tip}':/('${tip}'${var}:/" $outfile


done
echo End\; >> $outfile 
sed -i 's/index/\&index/g' $outfile
