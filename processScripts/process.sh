tips=$1
((tip_1=tips-1))
size=$2

rm -f true_anc_state
for ((i =1; i <=size; i++)) 
do
	awk -F , '{printf("%s,", $2)}' simulate/out.${i}.anc_state.csv  >> true_anc_state

	echo >> true_anc_state

done
awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var; print}' true_anc_state > true_anc_state.csv
rm true_anc_state.csv
