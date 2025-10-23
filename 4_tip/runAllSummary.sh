size=2500
tips=4

cd 3_bin
../fixOrder.sh
./sum.sh

cd ../8_cat
../fixOrder.sh

cd ../

cd 3_bin
/storage/anna/asr_phyddle/processScripts/process.sh $tips $size
/storage/anna/asr_phyddle/processScripts/parseParallel.sh $tips $size
#/storage/anna/asr_phyddle/mlParseParallel.sh $tips $size 
cd ../
