idx=1

nstate=8

outfile=empirical/lio.${idx}.ann_marg_daugh.nex
infile=empirical/lio.${idx}.form.tre
nodeLabel=empirical/lio.${idx}.node_labels.csv
dat=empirical/lio.${idx}.dat.csv
est=estimate/sim.empirical_est.labels_cat.csv

Rscript print_daughter_names.R
sed -i "s/'//g" empirical/parent_daughter.csv

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
	echo original name $nodeOld
	echo new name $nodeNew
	
	stateString="index=${index}"
	probString=
	
	i=1
	# Is this off by one?? No
	((col=(nodeNew)*nstate+2))

	# Order the probabilities, print the state number next to the probabiliites 
	probAll=$(grep ^${idx}, $est | awk -F , -v var1=$col -v var2=$nstate '{ for (i = 0; i < var2; i++) printf("%s\n", $(i+var1))}' )
	
	probA_B=$(echo $probAll | cut -d ' ' -f 1-2)
	probAB=$(echo $probAll | cut -d ' ' -f 3-8 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')

	ordered=$(echo $probA_B $probAB |  awk  '{ for (i = 0; i < 3; i++) printf("%s,%d\n", $(i+1), i)}' | sort -g -r)

	# Combines annoations for one node
	for line2 in $(echo $ordered |cut -d ' ' -f 1-3) 
	do
		
		prob=$(echo $line2 | awk -F ,  '{print $1}')
		state=$(echo $line2 | awk -F ,  '{print $2}')
		stateString="${stateString},end_state_${i}=${state}"
		probString="${probString},end_state_${i}_pp=${prob}"
		((i=i+1))
		
	done

	#other=$(echo $ordered | sed 's/ /\n/g' | tail -n+4 | awk -F',' '{sum+=$1} END{print sum;}')
	other=0


	## This will only work for up to three states
	probString="${probString},end_state_other_pp=${other}"


	# Find the line the original number is a daughter
	line3=$(cut -f 2-3 -d , empirical/parent_daughter.csv | grep -n  -e ${nodeOld}, -e ${nodeOld}$ | cut -f 1 -d :)

	# Find the parent of the original node number
	parent=$(sed "${line3}q;d" empirical/parent_daughter.csv | cut -f 1 -d , )
	
	# Find the parent with the relabeling 
	nodeParentNew=$(grep ${parent}, $nodeLabel | cut -f 2 -d , )

	echo original parent $parent
	echo new parent $nodeParentNew

	# ANNA_ check column -maybe wrong? 
	((col=(nodeParentNew)*nstate+2))
	head -1 $est | awk -F , -v var1=$col -v var2=$nstate '{ for (i = 0; i < var2; i++) printf("%s\n", $(i+var1))}' 
	
	# Find the column of the parent in the estimation file
	probAll=$(grep ^${idx}, $est | awk -F , -v var1=$col -v var2=$nstate '{ for (i = 0; i < var2; i++) printf("%s\n", $(i+var1))}' )
	echo $probAll

	left=$(cut -f 2-3 -d , empirical/parent_daughter.csv | grep  ${nodeOld} | cut -f 1 -d ,  )
	

	# If the tip of interests was the left daughter
	if [[ "$left" == "$nodeOld" ]]
	then
		#echo left

		start_probA=$(echo $probAll | cut -d ' ' -f 1,3,6 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probB=$(echo $probAll | cut -d ' ' -f 2,4,8 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probAB=$(echo $probAll | cut -d ' ' -f 5,7 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')

	else 
		#echo right 
		start_probA=$(echo $probAll | cut -d ' ' -f 1,4,5 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probB=$(echo $probAll | cut -d ' ' -f 2,3,7 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probAB=$(echo $probAll | cut -d ' ' -f 6,8 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
	fi

	echo $start_probA $start_probB $start_probAB

	ordered=$(echo $start_probA $start_probB  $start_probAB |  awk  '{ for (i = 0; i < 3; i++) printf("%s,%d\n", $(i+1), i)}' | sort -g -r)
	echo $ordered

	i=1
	d_stateString=
	d_probString=
	# Combines annoations for one node
	for line4 in $(echo $ordered |cut -d ' ' -f 1-3) 
	do
		
		prob=$(echo $line4 | awk -F ,  '{print $1}')
		state=$(echo $line4 | awk -F ,  '{print $2}')
		d_stateString="${d_stateString},start_state_${i}=${state}"
		d_probString="${d_probString},start_state_${i}_pp=${prob}"
		((i=i+1))
		
	done

	echo $d_stateString $d_probString
	# This will only work for up to three states
	d_probString="${d_probString},start_state_other_pp=${other}"

	sed -i -e "s/)${nodeOld}:/)[${stateString}${probString}${d_stateString}${d_probString}]:/" -e "s/)${nodeOld};/)[${stateString}${probString}${d_stateString}${d_probString}];/" $outfile

	
	((index=index+1))
#	exit
done

# Create tip annotations
# This will depend on the format of the input data
index=1
for line in $(cat ${dat} | tail -n+2) 
do
	echo $line
	tip=$(echo $line| awk -F , '{print $1}')
	echo $tip
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

	# Find the line the original number is a daughter
	line2=$(cut -f 2-3 -d , empirical/parent_daughter.csv | grep -n  -e ${tip}, -e ${tip}$ | cut -f 1 -d :)
	echo original num $line2

	# Find the parent of the original node number
	parent=$(sed "${line2}q;d" empirical/parent_daughter.csv | cut -f 1 -d , )
	echo parent $parent
	
	# Find the parent with the relabeling 
	nodeParentNew=$(grep $parent $nodeLabel | cut -f 2 -d , )
	echo parent new $nodeParentNew

	((col=(nodeParentNew)*nstate+2))
	
	# Find the column of the parent in the estimation file
	probAll=$(grep ^${idx}, $est | awk -F , -v var1=$col -v var2=$nstate '{ for (i = 0; i < var2; i++) printf("%s\n", $(i+var1))}' )

	left=$(cut -f 2-3 -d , empirical/parent_daughter.csv | grep  ${tip} | cut -f 1 -d ,  )
	

	# If the tip of interests was the left daughter
	if [[ "$left" == "$tip" ]]
	then
		echo left

		start_probA=$(echo $probAll | cut -d ' ' -f 1,3,6 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probB=$(echo $probAll | cut -d ' ' -f 2,4,8 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probAB=$(echo $probAll | cut -d ' ' -f 5,7 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')

	else 
		echo right 
		start_probA=$(echo $probAll | cut -d ' ' -f 1,4,5 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probB=$(echo $probAll | cut -d ' ' -f 2,3,7 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
		start_probAB=$(echo $probAll | cut -d ' ' -f 6,8 | sed 's/ /\n/g' | awk  '{sum+=$1} END{print sum;}')
	fi


	ordered=$(echo $start_probA $start_probB  $start_probAB |  awk  '{ for (i = 0; i < 3; i++) printf("%s,%d\n", $(i+1), i)}' | sort -g -r)

	i=1
	d_stateString=
	d_probString=
	# Combines annoations for one node
	for line3 in $(echo $ordered |cut -d ' ' -f 1-3) 
	do
		
		prob=$(echo $line3 | awk -F ,  '{print $1}')
		state=$(echo $line3 | awk -F ,  '{print $2}')
		d_stateString="${d_stateString},start_state_${i}=${state}"
		d_probString="${d_probString},start_state_${i}_pp=${prob}"
		((i=i+1))
		
	done

	# This will only work for up to three states
	d_probString="${d_probString},start_state_other_pp=${other}"

	#ANNA you need to add parent state of tips
	
	var="[index=${index}${d_stateString}${d_probString},end_state_1=${state_1},end_state_2=${state_2},end_state_3=${state_3},end_state_1_pp=${prob1},end_state_2_pp=${prob2},end_state_3_pp=${prob3},end_state_other_pp=0.0]"
	((index=index+1))

	# Combine these
	sed -i -e "s/,'${tip}':/,'${tip}'${var}:/" -e "s/('${tip}':/('${tip}'${var}:/" $outfile


done
echo End\; >> $outfile 
sed -i 's/index/\&index/g' $outfile
