
rm empirical/*
rm -f truth.csv
for i in {1..2500}
do

	for j in {1..3}
	do
		((k=(i-1)*3+j))
		# This is the tree
		ln simulate/out.${i}.tre empirical/out.${k}.tre 

		# This is the tip states
		ln simulate/out.${i}.dat.csv empirical/out.${k}.dat.csv

		# This is the one you are already estimating
		echo  asr_node_label > empirical/out.${k}.labels.csv
		echo $j >> empirical/out.${k}.labels.csv
		#head -$j simulate/out.${i}.anc_state.csv | tail -n 1 >> empirical/out.${k}.labels.csv

		head -${j} simulate/out.${i}.anc_state.csv | tail -n 1 | awk -F , '{print $2}' >> truth.csv

	done
done
