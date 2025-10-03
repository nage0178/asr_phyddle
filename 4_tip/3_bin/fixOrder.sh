echo idx,original,new > order.csv
for i in {1..2500}
do
	# Prints idx,original_node,new_node 
	tail -n+2  simulate/out.${i}.node_labels.csv |  sed 's/node//g' | awk -v var=$i '{print var "," $0}' >> order.csv
done
