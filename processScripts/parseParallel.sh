tips=$1
((tip_1=tips-1))
size=$2

mkdir -p parse
echo | awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var;  for (i = 1; i <= NF; ++i) printf("asr_%s,", i);}'  >  new_test_truth.csv
sed -i 's/^/idx,/g' new_test_truth.csv
echo >> new_test_truth.csv

echo | awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var;  for (i = 1; i <= NF; ++i) printf("asr_0_%s,asr_1_%s,", i, i);}'  >  new_test_est.csv
sed -i 's/^/idx,/g' new_test_est.csv
echo >> new_test_est.csv

head -1 estimate/out.test_true.labels_cat.csv  | sed 's/,/\n/g' > headerTrueLab
head -1 estimate/out.test_est.labels_cat.csv | sed 's/,/\n/g' > headerEstLab

length=$(wc -l estimate/out.test_true.labels_cat.csv | awk '{print $1}')
echo $length
for ((line2=2; line2<=length; line2++))
do
	echo $line2
	/storage/anna/asr_phyddle/processScripts/parseOne.sh $line2 &
	if (( line2 % 100 == 0 ))
	then
		echo waiting $line2
		wait
	fi
done
wait

for ((idx=1; idx <=size; idx++)) 
do
	cat parse/${idx}test_truth.csv >> new_test_truth.csv
	cat parse/${idx}test_est.csv >> new_test_est.csv
done

((tip2_1=tips*2-1))
awk -v var=$tips 'BEGIN{FS=OFS=","} {NF=var; print}' new_test_truth.csv > truth.csv
awk -v var=$tip2_1 'BEGIN{FS=OFS=","} {NF=var; print}' new_test_est.csv > est.csv

rm header*Lab
