if [ ! -f out.1.org_labels.csv ] ; then
	mv out.$1.labels.csv out.$1.org_labels.csv
fi

Rscript parse_json.R out.$1.org_labels.csv out.$1.labels.csv out.$1.json
