all: group_alleles calc_distances

group_alleles:
	python ../group_common_alleles.py > common_allele_pairs.txt
calc_distances:
	python ../main.py --whole_protein --batch common_allele_pairs.txt > common_distances_whole_protein.tsv
	python ../main.py --batch common_allele_pairs.txt > common_distances_binding_grooves.tsv

distances_somaticandtreatment:
	python real_hla_genotypes.py
	python ../main.py --batch real_allele_pairs.txt > real_distances_binding_grooves.tsv

#Make plots for HLA-B using common alleles and alleles occuring in somatic
plot:
	python plot.py --distances common_distances_binding_grooves.tsv --prefix common_binding_groove_ B
	python plot.py --distances real_distances_binding_grooves.tsv --prefix real_binding_groove_ B