#echo idx,asr_old_1,asr_old_2,asr_old_3,asr_old_4,asr_old_5,asr_old_6,asr_old_7,asr_old_8,asr_old_9 > new_test_truth.csv
#echo idx,asr_old_1_0,asr_old_1_1,asr_old_2_0,asr_old_2_1,asr_old_3_0,asr_old_3_1,asr_old_4_0,asr_old_4_1,asr_old_5_0,asr_old_5_1,asr_old_6_0,asr_old_6_1,asr_old_7_0,asr_old_7_1,asr_old_8_0,asr_old_8_1,asr_old_9_0,asr_old_9_1 > new_test_est.csv


echo | awk 'BEGIN{FS=OFS=","} {NF=99;  for (i = 1; i <= NF; ++i) printf("asr_%s,", i);}'  >  new_test_truth.csv
sed -i 's/^/idx,/g' new_test_truth.csv
echo >> new_test_truth.csv

echo | awk 'BEGIN{FS=OFS=","} {NF=99;  for (i = 1; i <= NF; ++i) printf("asr_0_%s,asr_1_%s,", i, i);}'  >  new_test_est.csv
sed -i 's/^/idx,/g' new_test_est.csv
echo >> new_test_est.csv

head -1 estimate/out.test_true.labels_cat.csv  | sed 's/,/\n/g' > headerTrueLab
head -1 estimate/out.test_est.labels_cat.csv | sed 's/,/\n/g' > headerEstLab

for line2 in $(tail +2 estimate/out.test_true.labels_cat.csv  )
#for line2 in $(tail +2 estimate/out.test_true.labels_cat.csv )
do
	idx=$(echo $line2 | awk -F , '{print $1}')
	echo idx  $idx

	oldLine=$idx
	newEstLine=$idx

	# This relies on this ordering matching the label ordering in the file
	for line in $( tail +2 simulate/out.${idx}.node_labels.csv | sed 's/node//g' | sed 's/,/ /g' | sort -k 1  -n | sed 's/^/node/g' | sed 's/ /,/g')
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
	echo $oldLine >> new_test_truth.csv
	echo $newEstLine >> new_test_est.csv
done

awk 'BEGIN{FS=OFS=","} {NF=100; print}' new_test_truth.csv > truth.csv
awk 'BEGIN{FS=OFS=","} {NF=199; print}' new_test_est.csv > est.csv
