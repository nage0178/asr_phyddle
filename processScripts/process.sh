#echo tipA,tipB,tipC,tipD,tipE,tipF,tipG,tipH,tipI,tipJ > tip_state
tips=$1
((tip_1=tips-1))
size=$2

rm -f  ml_asr
rm -f  ml_asr_prob
rm -f true_anc_state
for ((i =1; i <=size; i++)) 
do
	awk -F , '{printf("%s,", $2)}' simulate/out.${i}.asr.csv  >> ml_asr
	awk -F , '{printf("%s,", $4)}' simulate/out.${i}.asr.csv  >> ml_asr_prob
	awk -F , '{printf("%s,", $2)}' simulate/out.${i}.anc_state.csv  >> true_anc_state

	echo >> ml_asr
	echo >> ml_asr_prob
	echo >> true_anc_state

done
awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var; print}' ml_asr > ml_asr.csv
awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var; print}' ml_asr_prob > ml_asr_prob.csv
awk -v var=$tip_1 'BEGIN{FS=OFS=","} {NF=var; print}' true_anc_state > true_anc_state.csv
rm ml_asr
#rm ml_asr_prob
rm true_anc_state 
