sed -e 's/Lio/\nLio/g' -e 's/\[/\n\[/g'  data/small_tree.tre | grep Lio > small_taxa
grep -f small_taxa data/ranges.nex | sed 's/^ *//g' | sed -e 's/ \+/,/' -e 's/1$/,1/' -e 's/0$/,0/' > phyddle/empirical/lio.1.dat.csv
sed -i '1s/^/taxa,region1,region2\n/' phyddle/empirical/lio.1.dat.csv
grep Lio   data/small_tree.tre > data/lio.tre
#sed 's/\[\([^]]*\)\]//g'   data/small_tree.tre > data/lio_no_HPD.tre
Rscript add_node_annot.R

grep Lio data/ranges.nex | grep -v -f small_taxa > taxa_rm

grep -v -f taxa_rm data/ranges.nex > data/small_ranges.nex

sed -i 's/NTAX=199/NTAX=52/g' data/small_ranges.nex
