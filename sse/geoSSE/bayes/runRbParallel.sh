rb=~/revbayes/projects/cmake/build/rb
idx_max=2500
nthread=2500
mkdir -p screen_out

if ((idx_max % nthread == 0)) 
then
	((reps=idx_max/nthread))
else
	((reps=idx_max/nthread + 1))
fi
#((reps=idx_max/nthread + 1))

for ((i=1;i<nthread+1;i++))
do
	./runSetRb.sh $i $reps $idx_max &> screen_out/thread_$i & 

done
wait

