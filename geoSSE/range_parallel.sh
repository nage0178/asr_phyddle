# e.g. (1,0,0) means a species present only in region A, so it has range {A}.
# ranges_inv = {
#     (1,0,0):1,
#     (0,1,0):2,
#     (0,0,1):3,
#     (1,1,0):4,
#     (1,0,1):5,
#     (0,1,1):6,
#     (1,1,1):7
# }

#for i in {1..50000}
cd simulate
for i in {1..50000}
#for i in {6..7}
do
	../../one_range.sh $i & 
	if (( i % 100 == 0 ))
	then
		echo waiting $i
		wait
	fi
done
