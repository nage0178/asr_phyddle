mkdir -p data
#for idx in {1..10} 
for idx in {1..2500} 
do
	echo $idx
	in=../simulate/sim.${idx}.dat.csv
	out=data/out.${idx}.dat.nex
	echo \#NEXUS > $out
	echo Begin data\; >> $out
	echo Dimensions ntax=50 nchar=2\; >> $out
	echo Format datatype=Standard symbols=\"01\"\; >>  $out
	echo Matrix >>  $out
	
	tail -n+2 ${in} |  sed 's/,/ /'  | sed 's/,//' >> $out
	echo \; >> $out
	echo End\; >> $out
done


#out=data/test.${idx}.nex
#tree=data/test.${idx}.tre
##echo \#NEXUS > $out
#
#echo BEGIN TAXA\; > $out
##echo BEGIN TAXA\; >> $out
#echo        DIMENSIONS NTAX = 4\; >> $out
#echo        TAXLABELS >> $out
#
#tail -n+2 ${in} | cut -d , -f 1  >> $out
#echo	\; >> $out
#echo END\; >> $out
#echo BEGIN TREES\; >> $out
#cat $tree >> $out
#echo END\; >> $out



