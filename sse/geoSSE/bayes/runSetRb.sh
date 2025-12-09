RB_CMD=~/revbayes/projects/cmake/rb

thread=$1
reps=$2
max=$3

for ((i=1;i<=reps;i++))
do

	if [[ $max -lt $idx ]] ; then
		exit
	fi

	((idx=(thread-1)*reps+i))
	if [ ! -f output/1/${idx}_ase.tre ] ; then
	#	echo $idx found
		RB_STR="IDX=${idx}; run_i=1; source(\"scripts/geosse.Rev\")"
		#echo ${RB_STR} $RB_CMD
		#echo ${RB_STR} | ${RB_CMD} &> screen_out/${idx}_1.txt & 
		echo ${RB_STR} | ${RB_CMD} &> screen_out/${idx}_1.txt 
		echo $idx 1
	fi


	#$rb scripts/mcmc_BiSSE_root_stat.Rev ${idx} &> screen_out/${idx}.txt 
	
	

	if [ ! -f output/2/${idx}_ase.tre ] ; then
		RB_STR="IDX=${idx}; run_i=2; source(\"scripts/geosse.Rev\")"
		#echo ${RB_STR} $RB_CMD
		#echo ${RB_STR} | ${RB_CMD} &> screen_out/${idx}_2.txt  & 
		echo ${RB_STR} | ${RB_CMD} &> screen_out/${idx}_2.txt  
		echo $idx 2
	fi

done

wait
