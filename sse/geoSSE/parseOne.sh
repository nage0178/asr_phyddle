lineNum=$1
pre=sim
line2=$(head -$lineNum estimate/${pre}.test_true.labels_cat.csv | tail -n 1)
idx=$(echo $line2 | awk -F , '{print $1}')


oldLine=$idx
newEstLine=$idx

# This relies on this ordering matching the label ordering in the file
for line in $( tail +2 simulate/${pre}.${idx}.node_labels.csv | sed 's/nd//g' | sed 's/,/ /g' | sort -k 1  -n | sed 's/^/nd/g' | sed 's/ /,/g')
do
	original=$(echo $line |awk -F ,  '{print $1}')
	new=$(echo $line | awk -F ,  '{print $2}')


	# Find column with matching name in true
	colTrue=$( grep -n asr_$new$ headerTrueLab | awk -F : '{print $1}')

	# find the true value in the original node labels file, not the inference output
	trueVal=$(echo $line2 | awk -F , -v num=$colTrue '{print $num}')
	
	oldLine=${oldLine},${trueVal}

	# Find column with matching name in test
	colEst=$( grep -n asr_${new}_ headerEstLab | awk -F : '{print $1}')

	cols=$(echo $colEst)
	est=$(grep ^${idx}, estimate/${pre}.test_est.labels_cat.csv | cut -d , -f"${cols}"  )

	newEstLine=${newEstLine},${est}


done
echo $oldLine > parse/${idx}test_truth.csv
echo $newEstLine > parse/${idx}test_est.csv
