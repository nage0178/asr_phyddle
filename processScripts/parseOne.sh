lineNum=$1
line2=$(head -$lineNum estimate/out.test_true.labels_cat.csv | tail -n 1)

idx=$(echo $line2 | awk -F , '{print $1}')
#echo idx  $idx


oldLine=$idx
newEstLine=$idx

# This relies on this ordering matching the label ordering in the file
for line in $( tail +2 simulate/out.${idx}.node_labels.csv | head -2500 | sed 's/node//g' | sed 's/,/ /g' | sort -k 1  -n | sed 's/^/node/g' | sed 's/ /,/g')
do
	original=$(echo $line |awk -F ,  '{print $1}')
	new=$(echo $line | awk -F ,  '{print $2}')


	# Find column with matching name in true
	colTrue=$( grep -n asr_$new$ headerTrueLab | awk -F : '{print $1}')
	#colTrue=$(head -1 estimate/out.test_true.labels_cat.csv  | sed 's/,/\n/g'  | grep -n asr_$new$ | awk -F : '{print $1}')

	# Here we are finding the true value in the original node labels file, not the inference output
	trueVal=$(echo $line2 | awk -F , -v num=$colTrue '{print $num}')
	
	oldLine=${oldLine},${trueVal}

	# Find column with matching name in test
	colEst=$( grep -n asr_${new}_ headerEstLab | awk -F : '{print $1}')
	#colEst=$(head -1 estimate/out.test_est.labels_cat.csv | sed 's/,/\n/g' | grep -n asr_${new}_ | awk -F : '{print $1}')
	colZero=$(echo $colEst | awk  '{print $1}')
	colOne=$(echo $colEst | awk  '{print $2}')

	estZero=$(grep ^${idx}, estimate/out.test_est.labels_cat.csv | awk  -v var="$colZero" -F , '{print $var }' )
	estOne=$(grep ^${idx}, estimate/out.test_est.labels_cat.csv | awk  -v var="$colOne" -F , '{print $var }' )

	newEstLine=${newEstLine},${estZero},${estOne}


done
echo $oldLine > parse/${idx}test_truth.csv
echo $newEstLine > parse/${idx}test_est.csv
