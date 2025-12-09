rb=~/revbayes/projects/cmake/build/rb
thread=$1
reps=$2
max=$3

for ((i=1;i<=reps;i++))
do
	((idx=(thread-1)*reps+i))
	if [ -f output/${idx}_BiSSE_stoch_map_character.tree ] ; then
		echo $idx found
		continue
	fi

	if [[ $max -lt $idx ]] ; then
		exit
	fi
	echo $idx
	$rb scripts/mcmc_BiSSE_root_stat.Rev ${idx} &> screen_out/${idx}.txt 

done

wait
