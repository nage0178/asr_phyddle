size=$1
echo idx,rate > rates.csv

for ((i =1; i <=size; i++))
do
	rate=$(cat simulate/out.${i}.rate.csv)
	echo $i,$rate  >> rates.csv
done
