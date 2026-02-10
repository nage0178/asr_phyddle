tips=50
size=2500
dir=../../processScripts/

${dir}findTimes.sh $tips $size
${dir}process.sh $tips $size
${dir}parseParallel.sh $tips $size
