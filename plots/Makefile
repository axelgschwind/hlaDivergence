all: group_alleles calc_distances

group_alleles:
	python ../group_common_alleles.py > allele_pairs.txt
calc_distances:
	python ../main.py --batch allele_pairs.txt > distances.tsv