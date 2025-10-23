#rm -r accuracy
phyddle -s SFTEP &> outAll
mv train train_first_train
mv estimate estimate_first_train

# Simulate dataset
phyddle -s S --sim_dir accuracy -c config_empirical.py

#Create files with labels reversed
for i in {150001..175000}
do
	((j=i+25000))
	sed 's/,0/,2/g' accuracy/out.${i}.anc_state.csv > accuracy/out.${j}.anc_state.csv
	sed -i 's/,1/,0/g' accuracy/out.${j}.anc_state.csv
	sed -i 's/,2/,1/g' accuracy/out.${j}.anc_state.csv

	sed 's/,0/,2/g' accuracy/out.${i}.dat.csv > accuracy/out.${j}.dat.csv
	sed -i 's/,1/,0/g' accuracy/out.${j}.dat.csv
	sed -i 's/,2/,1/g' accuracy/out.${j}.dat.csv
	

	sed 's/,0/,2/g' accuracy/out.${i}.labels.csv > accuracy/out.${j}.labels.csv
	sed -i 's/,1/,0/g' accuracy/out.${j}.labels.csv
	sed -i 's/,2/,1/g' accuracy/out.${j}.labels.csv

	cp accuracy/out.${i}.tre accuracy/out.${j}.tre 
done

# Format and estimate with phyddle
phyddle -s F --no_sim --emp_dir accuracy -c config_empirical.py
phyddle -s E --no_sim --emp_dir accuracy -c config_empirical.py

## Format 
#mkdir -p format_accuracy/simulate
#cd format_accuracy
#phyddle -s F --sim_dir ../accuracy --prop_test 0  -c ../config_empirical.py
#phyddle -s F --sim_dir ../accuracy --prop_test 0  -c ../config_empirical.py
