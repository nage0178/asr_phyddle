mkdir -p data
for idx in {1..2500} 
do
	ntax=$(tail -n+2 ../simulate/out.${idx}.dat.csv | wc -l)
	echo $idx
	in=../simulate/out.${idx}.dat.csv
	out=data/out.${idx}.dat.nex
	echo \#NEXUS > $out
	echo Begin data\; >> $out
	echo Dimensions ntax=${ntax} nchar=1\; >> $out
	echo Format datatype=Standard symbols=\"01\"\; >>  $out
	echo Matrix >>  $out
	
	tail -n+2 ${in} |  sed 's/,/ /g' >> $out
	echo \; >> $out
	echo End\; >> $out
done
