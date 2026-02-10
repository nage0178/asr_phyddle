size=2500
dir=../../../processScripts/
#dir=/storage/anna/asr_phyddle/processScripts/
tips=50

${dir}findTimes.sh $tips $size
${dir}process.sh $tips $size
${dir}parseParallel.sh $tips $size
${dir}mlParseParallel.sh $tips $size 
