for dir in fix var
do

	for size in 50 100 200 
	do
		mkdir -p $dir/$size


		if [ "$dir" =  "fix" ];
		then
		 	sed "s/nupper/$size/g" config.py > $dir/$size/config.py
			sed -i "s/nlower/$size/g"  $dir/$size/config.py

		else
		 	sed "s/nupper/$size/g" config.py > $dir/$size/config.py
			((lower = size / 5))
			echo $lower
			sed -i "s/nlower/$lower/g"  $dir/$size/config.py
		fi
		
	done
done
