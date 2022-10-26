import itertools
import sys
from collections import defaultdict
import os

def main(argv):
	this_dir, this_filename = os.path.split(__file__)

	data_dir = os.path.join(this_dir, "data", "common_alleles.txt")

	common_alleles = defaultdict(list)
	with open( data_dir ) as file:
		for line in file:
			(locus, allele) = line.strip().split('*')
			common_alleles[locus].append(line.strip())

	for key in common_alleles.keys():
		for pair in itertools.combinations(common_alleles[key], 2):
			print(pair[0], pair[1], sep="\t")

if __name__ == "__main__":
	main(sys.argv)