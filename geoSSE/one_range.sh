i=$1
awk -F , '{print $1 "," $4 }' cladogenetic\ event\ history_${i}.csv > range_${i}.csv
#pwd
#ls cladogenetic\ event\ history_${i}.csv 

sed -i 's/idx,final_state/idx,final_state,region_1,region_2,region_3/' range_$i.csv
sed -i 's/1$/1,1,0,0/g' range_${i}.csv
sed -i 's/2$/2,0,1,0/g' range_${i}.csv
sed -i 's/3$/3,0,0,1/g' range_${i}.csv
sed -i 's/4$/4,1,1,0/g' range_${i}.csv
sed -i 's/5$/5,1,0,1/g' range_${i}.csv
sed -i 's/6$/6,0,1,1/g' range_${i}.csv
sed -i 's/7$/7,1,1,1/g' range_${i}.csv

tail -n+2  range_${i}.csv | awk -F , '{print "node" $1 "," $3}'   > sim.${i}.anc_state.csv
