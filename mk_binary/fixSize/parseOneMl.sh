lineNum=$1
line2=$(head -$lineNum estimate/out.test_true.labels_cat.csv | tail -n 1)

idx=$(echo $line2 | awk -F , '{print $1}')
#echo idx  $idx


oldLine=$idx
newEstLine=$idx

# This relies on this ordering matching the label ordering in the file
for line in $( tail +2 simulate/out.${idx}.node_labels.csv)
do
	original=$(echo $line |awk -F ,  '{print $1}')
	new=$(echo $line | awk -F ,  '{print $2}')


	# Find column with matching name in true
#	colTrue=$( grep -n asr_$new headerTrueLab | awk -F : '{print $1}')
	estML=$(grep ${original}, simulate/out.${idx}.asr.csv | awk -F , '{print $2}')


	newEstLine=${newEstLine},${estML}


done
echo $newEstLine > parse/${idx}test_estML.csv
