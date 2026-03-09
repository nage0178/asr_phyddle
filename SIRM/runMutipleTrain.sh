phyddle -s S &> outSim 
phyddle -s T  &> outT1
phyddle -s T --trn_dir train_2 --est_dir estimate_2 --plt_dir plot_2 &> outT2
phyddle -s T --trn_dir train_3 --est_dir estimate_3 --plt_dir plot_3 &> outT3 
phyddle -s T --trn_dir train_4 --est_dir estimate_4 --plt_dir plot_4 &> outT4 
phyddle -s T --trn_dir train_5 --est_dir estimate_5 --plt_dir plot_5 &> outT5 


# These are using prevalence (2 weeks)
phyddle -s E --trn_dir train --est_dir estimate --plt_dir plot --no_sim &> outE1
cp empirical/out.1.est.tre empirical/train_1_out.1.est.tre

phyddle -s E --trn_dir train_2 --est_dir estimate_2 --plt_dir plot_2 --no_sim &> outE2
cp empirical/out.1.est.tre empirical/train_2_out.1.est.tre

phyddle -s E --trn_dir train_3 --est_dir estimate_3 --plt_dir plot_3 --no_sim &> outE3 
cp empirical/out.1.est.tre empirical/train_3_out.1.est.tre

phyddle -s E --trn_dir train_4 --est_dir estimate_4 --plt_dir plot_4 --no_sim &> outE4 
cp empirical/out.1.est.tre empirical/train_4_out.1.est.tre

phyddle -s E --trn_dir train_5 --est_dir estimate_5 --plt_dir plot_5 --no_sim &> outE5 
cp empirical/out.1.est.tre empirical/train_5_out.1.est.tre

