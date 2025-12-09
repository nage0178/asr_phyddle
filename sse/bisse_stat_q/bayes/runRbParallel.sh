idx_max=2500
nthread=75
((reps=idx_max/nthread + 1))
for ((i=1;i<nthread+1;i++))
do
	./runSetRb.sh $i $reps $idx_max &> screen_out/thread_$i & 

done
wait

