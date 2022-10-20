import sys
from optparse import OptionParser
from typing import Sequence
import pandas as pd
from Bio import AlignIO, SeqRecord
from Bio.Seq import Seq


#Brings allele id to the form A*01:01, i.e. nomenclature for 4 digits without "HLA-"
def parse_allele_name(allele: str):
	(locus, parts) = allele.replace("HLA-","").split("*")
	parts = parts.split(":")
	if len(parts) == 0 or len(parts) > 4:
		raise Exception("Cannot parse HLA string " + allele)
	if len(parts) == 1:
		parts.append("01")

	return locus + "*" + ":".join(parts[0:2])


#exon coordinates are compatible with IMGT/HLA protein aligments. Exon 2/3 form binding grooves of MHC class I. Exon 2 forms binding groove of MHC class II
_binding_groove_coords = {
	"A": slice(26,228),
	"B": slice(30,241),
	"C": slice(24,235),
	"DQA1": slice(27,111),
	"DQB1": slice(36,127),
	"DRB1": slice(33,127)
}
#return first allele that matches four digits allele
def get_protein_sequence(allele, alignments: AlignIO.MultipleSeqAlignment, whole_protein: bool=False):
	allele = parse_allele_name(allele)
	for alignment in alignments:
		if allele == parse_allele_name(alignment.name):
			if whole_protein:
				return alignment.seq
			return alignment.seq[ _binding_groove_coords[get_hla_locus(allele)] ]

	raise Exception("Could not find alignment for %s" %allele)


#Calculate average Grantham distance for two sequences of equal length
_grantham_matrix = pd.DataFrame()
def grantham_distance(seq1: Seq, seq2: Seq):
	#initialize grantham matrix with values stored in grantham_matrix.tsv
	global _grantham_matrix
	if len(_grantham_matrix) == 0:
		_grantham_matrix = pd.read_csv("data/grantham.tsv", sep="\t", header=0, index_col=0, comment='#')

	if len(seq1) != len(seq2):
		raise Exception("Alignemnts of equal length are neccessary to calculate Grantham distance!")

	total_length = 0
	aa_distance = 0

	for index,aa1 in enumerate(seq1):
		aa2 = seq2[index]
		#ignore unknown amino acids (i.e. insertions and deletions)
		if aa1 not in _grantham_matrix.columns or aa2 not in _grantham_matrix.columns:
			continue

		total_length += 1

		if aa1 != aa2:
			aa_distance += _grantham_matrix[aa1][aa2]

	return aa_distance/total_length

def get_hla_locus(name: str):
	parts = name.replace("HLA-","").split("*")
	return parts[0]



def main(argv):
	parser = OptionParser(usage="usage: %prog [options] HLA-ALLELE1 HLA-ALLELE2", description="Calculates HLA diversity metrics. Pass allele names according HLA nomenclature.")
	parser.add_option("--whole_protein", action="store_true", default=False, dest="whole_protein", help="Use the whole protein sequence of MHC complex. Otherwise only binding grooves will be considered (Exon 2/3 for MHC class I and exon 2 for class II, respectively).")

	(options, args) = parser.parse_args()
	
	if(len(args) != 2):
		raise Exception("Specify two HLA alleles!")
	allele1 = args[0]
	allele2 = args[1]

	#aligned AA sequences

	locus = get_hla_locus(allele1)
	locus2 = get_hla_locus(allele2)
	if locus != locus2:
		raise Exception("The loci HLA-%s and HLA-%s must be the same" % (locus, locus2))

	#IMGT/HLA protein aligments
	protein_alignments = AlignIO.read("data/%s_prot.msf" % locus, "msf")

	seq1 = get_protein_sequence(allele1, protein_alignments, options.whole_protein)
	seq2 = get_protein_sequence(allele2, protein_alignments, options.whole_protein)
	print(grantham_distance(seq1, seq2))

if __name__ == "__main__":
	main(sys.argv)
