tips=$1
((tip_1=tips-1))
size=$2

mkdir -p parse

echo | awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var;  for (i = 1; i <= NF; ++i) printf("asr_%s,", i);}'  >  new_ml_est.csv
sed -i 's/^/idx,/g' new_ml_est.csv
echo >> new_ml_est.csv


length=$(wc -l estimate/out.test_true.labels_cat.csv | awk '{print $1}')
echo $length
for ((line2=2; line2<=length; line2++))
do
	echo $line2
	../parseOneMl.sh $line2 &
	if (( line2 % 100 == 0 ))
	then
		echo waiting $line2
		wait
	fi
done
wait

for ((idx=1; idx <=size; idx++)) 
do
	cat parse/${idx}test_estML.csv >> new_ml_est.csv
done

awk -v var=$tips 'BEGIN{FS=OFS=","} {NF=var; print}' new_ml_est.csv > ml_est_reorder.csv
