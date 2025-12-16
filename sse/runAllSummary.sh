tips=50
size=2500
dir=/storage/anna/asr_phyddle/processScripts/

cd bisse_stat_q
${dir}findTimes.sh $tips $size
${dir}process.sh $tips $size
${dir}parseParallel.sh $tips $size
cd ../
