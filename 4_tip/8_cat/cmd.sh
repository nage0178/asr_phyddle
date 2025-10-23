mkdir simulate
# Train and estimate with first network
phyddle -s FTEP --sim_dir ../3_bin/simulate &> outFTEP & 

# Run additional simulations as empirical data
phyddle -s F --no_sim --emp_dir ../3_bin/accuracy -c config_empirical.py
mkdir -p accuracy_est
phyddle -s E --no_sim --est_dir accuracy_est -c config_empirical.py

# Train second network
phyddle -s T --trn_prefix train2 &> outTrain 
phyddle -s E --trn_prefix train2 --est_prefix est2 &> outEst2 
