
# Labels internal nodes
# Also need to remove negative branch lengths and subsample
Rscript label_tree.R


echo \#NEXUS > out.1.dat.nex
echo Begin DATA\; >> out.1.dat.nex
echo Dimensions NTAX=50 NCHAR=5\; >> out.1.dat.nex
#echo Dimensions NTAX=100 NCHAR=5\; >> out.1.dat.nex
echo Format MISSING=? GAP=- DATATYPE=STANDARD SYMBOLS="01"\; >> out.1.dat.nex
echo Matrix >> out.1.dat.nex

# Will need to subsample with tips that are acutally used
grep "|" Ebola_glm_mascot.mcc.trees | grep -E  " [1-9][0-9]?[0-9]? " | awk -F "|" '{print $3 "," $5}' | sed 's/ EBOV//g' | sed -e 's/\t// g' -e 's/ //g'  | sed 's/,/\t/g' > locations.txt

#sed -i -e 's/Kailahun/0/g' -e 's/Kenema/0/g'  -e 's/Kono/1/g' -e 's/Koinadugu/1/g' -e 's/Bombali/1/g' -e 's/Tonkolili/2/g' -e 's/PortLoko/2/g' -e 's/Kambia/2/g' -e 's/Pujehun/3/g'  -e 's/Moyamba/3/g' -e 's/Bonthe/3/g' -e 's/WesternRural/4/g' -e 's/WesternUrban/4/g' locations.txt 
#sed -i 's/Bo/0/g' locations.txt

sed -i -e 's/Kailahun/10000/g' -e 's/Kenema/10000/g'  -e 's/Kono/01000/g' -e 's/Koinadugu/01000/g' -e 's/Bombali/01000/g' -e 's/Tonkolili/00100/g' -e 's/PortLoko/00100/g' -e 's/Kambia/00100/g' -e 's/Pujehun/00010/g'  -e 's/Moyamba/00010/g' -e 's/Bonthe/00010/g' -e 's/WesternRural/00001/g' -e 's/WesternUrban/00001/g' locations.txt 
sed -i 's/Bo/10000/g' locations.txt

grep  -f kept_samples.csv locations.txt >> out.1.dat.nex

echo \; >> out.1.dat.nex
echo END\; >> out.1.dat.nex
## Red= 0
#Kailahun
#Kenema
#Bo
#
##Purple = 1
#Kono
#Koinadugu
#Bombali
#
##Yellow = 2
#Tonkolili
#PortLoko
#Kambia
#
##Green = 3 
#Pujehun
#Moyamba
#Bonthe
#
## Blue = 4
#WesternRural
#WesternUrban
#
cp out.1.dat.nex out.1.dat.csv
