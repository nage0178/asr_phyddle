phyddle -s SFTEP &> outALL

# Estimate with larger dataset
mkdir -p accuracy/simulate/
phyddle -s F --emp_dir accuracy/empirical --no_sim -c config_empirical.py
phyddle -s E  --est_dir accuracy/estimate --no_sim -c config_empirical.py

# Second training
phyddle -s T --trn_prefix train2 &> outTrain & 
phyddle -s E --trn_prefix train2 --est_prefix est2 &> outEst2 & 
