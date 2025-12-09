rb=~/revbayes/projects/cmake/build/rb
mkdir -p screen_out
idx_max=2500
nthread=20
if ((idx_max % nthread == 0)) 
then
	((reps=idx_max/nthread))
else
	((reps=idx_max/nthread + 1))
fi

cd ${dir}/bayes
for ((i=1;i<nthread+1;i++))
	do
		./runSetRb.sh $i $reps $idx_max &> screen_out/thread_$i & 
	
	done
wait

