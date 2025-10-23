
#rm -r accuracy/
mkdir -p accuracy/empirical
for i in {150001..200000}
do

	for j in {1..3}
	do
		((k=(i-1)*3+j))
		# This is the tree
		ln ../3_bin/accuracy/out.${i}.tre accuracy/empirical/out.${k}.tre 

		# This is the tip states
		ln ../3_bin/accuracy/out.${i}.dat.csv accuracy/empirical/out.${k}.dat.csv

		# This is the one you are already estimating
		echo  asr_node_label > accuracy/empirical/out.${k}.labels.csv
		echo node$j >> accuracy/empirical/out.${k}.labels.csv


	done
done
