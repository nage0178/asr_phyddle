rb=~/revbayes/projects/cmake/build/rb
thread=$1
reps=$2
max=$3

for ((i=1;i<=reps;i++))
do
	((idx=(thread-1)*reps+i))
	if [ -f output/${idx}_ase.tree ] ; then
		echo $idx found
		continue
	fi

	if [[ $max -lt $idx ]] ; then
		exit
	fi
	echo $idx

	$rb scripts/mcmc_mk.Rev ${idx} &> screen_out/${idx}.txt 

done

wait
