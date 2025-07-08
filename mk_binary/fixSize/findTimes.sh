#echo tipA,tipB,tipC,tipD,tipE,tipF,tipG,tipH,tipI,tipJ > tip_state
tips=$1
((tip_1=tips-1))
size=$2

echo | awk -v var="$tip_1" 'BEGIN{FS=OFS=","} {NF=var;  for (i = 1; i <= NF; ++i) printf("node%s,", i);}'  > prop_ages

sed -i 's/node1,/idx,node1,/g' prop_ages
echo >> prop_ages

for ((i =1; i <=size; i++)) #{1..2500}
do
	#tail -n+2 simulate/out.${i}.dat.csv |awk -F , '{printf("%s,", $2)}' >> tip_state

	n_times=$(tail -n 1 simulate/out.${i}.age_prop.csv)
	echo ${i},${n_times}>> prop_ages 
	

done
awk -v var=$tips 'BEGIN{FS=OFS=","} {NF=var; print}' prop_ages > prop_ages.csv
