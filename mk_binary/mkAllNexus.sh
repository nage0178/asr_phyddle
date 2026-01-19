scrpt=../../../../processScripts
#scrpt=/storage/anna/asr_phyddle/processScripts
for dir in {fix,var}/{50,100,200} 
do
	cd $dir
	mkdir bayes
	cd bayes
	${scrpt}/mkNexusBayes.sh
	cd ../../../
	
done
