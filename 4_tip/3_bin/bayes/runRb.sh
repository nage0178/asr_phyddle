rb=~/revbayes/projects/cmake/build/rb

mkdir screen_out

for idx in {1..2500}
do
	$rb scripts/mcmc_4_tip.Rev ${idx} > screen_out/${idx}.txt &

	if (( idx % 100 == 0 ))
	then
		echo waiting $idx
		wait
	fi

done
