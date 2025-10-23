source ~/phyddle/devPhyddle/bin/activate

cd 50
for size in 50 
do
	rm -r empirical_${size}
	mkdir empirical_${size}
	
	cp ../../fix/${size}/simulate/out.{1..2500}.tre empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.dat.csv empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.anc_state.csv empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.labels.csv empirical_${size}/ & 
	wait
	
	phyddle -s F --no_sim --emp_dir empirical_${size}
	phyddle -s E --no_sim --est_prefix ${size}

done
cd ../

cd 100
for size in 50 100
do
	rm -r empirical_${size}
	mkdir empirical_${size}
	
	cp ../../fix/${size}/simulate/out.{1..2500}.tre empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.dat.csv empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.anc_state.csv empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.labels.csv empirical_${size}/ & 
	wait
	
	phyddle -s F --no_sim --emp_dir empirical_${size}
	phyddle -s E --no_sim --est_prefix ${size}

done
cd ../

cd 200 
for size in 50 100 200
do
	rm -r empirical_${size}
	mkdir empirical_${size}
	
	cp ../../fix/${size}/simulate/out.{1..2500}.tre empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.dat.csv empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.anc_state.csv empirical_${size}/ & 
	cp ../../fix/${size}/simulate/out.{1..2500}.labels.csv empirical_${size}/ & 
	wait
	
	phyddle -s F --no_sim --emp_dir empirical_${size}
	phyddle -s E --no_sim --est_prefix ${size}

done
cd ../
