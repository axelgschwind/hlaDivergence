all: group_alleles calc_distances plot

group_alleles:
	python ../group_common_alleles.py > allele_pairs.txt
calc_distances:
	python ../main.py --whole_protein --batch allele_pairs.txt > distances_whole_protein.tsv
	python ../main.py --batch allele_pairs.txt > distances_binding_grooves.tsv
plot:
	python plot.py