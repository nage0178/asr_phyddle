size=2500

for tips in 50 100 200
do 
	cd ${tips}
	/storage/anna/asr_phyddle/processScripts/findTimes.sh $tips $size
	/storage/anna/asr_phyddle/processScripts/process.sh $tips $size
        /storage/anna/asr_phyddle/processScripts/parseParallel.sh $tips $size
	/storage/anna/asr_phyddle/processScripts/mlParseParallel.sh $tips $size 
	cd ../
done
