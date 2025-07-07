mkdir empirical
cd empirical
for i in {1..2500}
do
	ln -s ../simulate/sim.${i}.tre .
	#ln -s simulate/sim.${i}.tre empirical/.

	ln -s ../simulate/sim.${i}.labels.csv .
	#ln -s simulate/sim.${i}.labels.csv empirical/.

	# Need to reorder the data
	echo taxa,region_1,region_2,region_3 > sim.${i}.dat.csv
	#echo taxa,region_1,region_2,region_3 > empirical/sim.${i}.dat.csv
	awk -F , '{print $1 "," $3 "," $2 "," $4 }' ../simulate/sim.${i}.dat.csv >> sim.${i}.dat.csv
	#awk -F , '{print $1 "," $3 "," $2 "," $4 }' simulate/sim.${i}.dat.csv >> empirical/sim.${i}.dat.csv

	
	tail -n+2  ../simulate/range_${i}.csv | awk -F , '{print "node" $1 "," $4}'   > sim.${i}.anc_state.csv
	#tail -n+2  simulate/range_${i}.csv | awk -F , '{print "node" $1 "," $4}'   > empirical/sim.${i}.anc_state.csv

done


