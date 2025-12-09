rb=~/revbayes/projects/cmake/build/rb
idx_max=2500
nthread=50
if ((idx_max % nthread == 0)) 
then
	((reps=idx_max/nthread))
else
	((reps=idx_max/nthread + 1))
fi

for dir in {fix,var}/{50,100,200}
do
	echo $dir
	cd ${dir}/bayes
	for ((i=1;i<nthread+1;i++))
		do
			../../../runSetRb.sh $i $reps $idx_max &> screen_out/thread_$i & 
		
		done
		cd ../../../
			wait
done
wait

