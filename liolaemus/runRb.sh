RB_CMD=~/revbayes/projects/cmake/rb


if [ ! -f output/1/${idx}_ase.tre ] ; then
	RB_STR="IDX=${idx}; run_i=1; source(\"scripts/geosse_match_phyddle.Rev\")"
	echo ${RB_STR} | ${RB_CMD} &> screen_out/lio_1.txt & 
	#echo ${RB_STR} | ${RB_CMD} &> screen_out/${idx}_1.txt 
fi

if [ ! -f output/2/${idx}_ase.tre ] ; then
	RB_STR="IDX=${idx}; run_i=2; source(\"scripts/geosse_match_phyddle.Rev\")"
	echo ${RB_STR} | ${RB_CMD} &> screen_out/lio_2.txt  & 
	#echo ${RB_STR} | ${RB_CMD} &> screen_out/${idx}_2.txt  
fi


wait
