size=2500
dir=/storage/anna/asr_phyddle/processScripts/
for tips in 50 100 200
do 
	cd ${tips}
	${dir}findTimes.sh $tips $size
	${dir}process.sh $tips $size
        ${dir}parseParallel.sh $tips $size
	${dir}mlParseParallel.sh $tips $size 
	cd ../
done
