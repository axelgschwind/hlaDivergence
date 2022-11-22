import sys
import os
from os.path import exists
from optparse import OptionParser
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq, MutableSeq
from scipy.spatial.distance import pdist, squareform
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord

#Brings allele id to the form A*01:01, i.e. nomenclature for 4 digits without "HLA-"
def parse_allele_name(allele: str):
	(locus, parts) = allele.strip().replace("HLA-","").split("*")
	parts = parts.split(":")
	if len(parts) == 0 or len(parts) > 5:
		raise Exception("Cannot parse HLA string " + allele)
	if len(parts) == 1:
		parts.append("01")

	return locus + "*" + ":".join(parts[0:2])

#exon coordinates are compatible with IMGT/HLA protein aligments. Exon 2/3 form binding grooves of MHC class I. Exon 2 forms binding groove of MHC class II
#exonic "borders" were taken from CCDS release 24 and exclude AA along splicing junctions. "Borders" are ajusted for the IMGT alignments
_binding_groove_coords = {
	"A": slice(27,227),
	"B": slice(31,239),
	"C": slice(25,232),
	"DQA1": slice(28,110),
	"DQB1": slice(37,126),
	"DRB1": slice(34,125)
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
		_grantham_matrix = pd.read_csv(data_path("grantham.tsv"), sep="\t", header=0, index_col=0, comment='#')

	if len(seq1) != len(seq2):
		raise Exception("Alignemnts of equal length are neccessary to calculate Grantham distance!")

	total_length = 0
	aa_distance = 0

	for index,aa1 in enumerate(str(seq1)):
		aa2 = seq2[index]
		#ignore unknown amino acids (i.e. insertions and deletions)
		if aa1 not in _grantham_matrix.columns or aa2 not in _grantham_matrix.columns:
			continue

		total_length += 1

		if aa1 != aa2:
			aa_distance += _grantham_matrix[str(aa1)][str(aa2)]

	return aa_distance/total_length

#Calculate average Sandberg distance for two sequences of equal length
_sandberg_matrix = pd.DataFrame()
def sandberg_distance(seq1: Seq, seq2: Seq):
	global _sandberg_matrix
	if(len(_sandberg_matrix) == 0):
		#get z values from Sandberg et al. original paper
		sandberg_z_values = pd.read_csv(data_path("sandberg.tsv"), sep="\t", header=0, index_col=0, comment='#')
		#calculate euclidian distance between the vectors of all amino acids
		_sandberg_matrix = pd.DataFrame(squareform( pdist(sandberg_z_values) ), index=sandberg_z_values.index.values, columns = sandberg_z_values.index.values )
	
	if len(seq1) != len(seq2):
		raise Exception("Alignemnts of equal length are neccessary to calculate Sandberg distance!")
	total_length = 0
	aa_distance = 0

	for index,aa1 in enumerate(str(seq1)):
		aa2 = seq2[index]
		#ignore unknown amino acids (i.e. insertions and deletions)
		if aa1 not in _sandberg_matrix.columns or aa2 not in _sandberg_matrix.columns:
			continue
		total_length += 1

		if aa1 != aa2:
			aa_distance += _sandberg_matrix[str(aa1)][str(aa2)]

	return aa_distance / total_length

def p_distance(seq1: Seq, seq2: Seq):
	global _grantham_matrix
	if(len(seq1) != len(seq2)):
		raise Exception("Alignemnts of equal length are neccessary to calculate pDistance!")

	total_length = 0
	aa_distance = 0

	for index,aa1 in enumerate(str(seq1)):
		aa2 = seq2[index]
		#ignore unknown amino acids (i.e. insertions and deletions)
		if aa1 not in _grantham_matrix.columns or aa2 not in _grantham_matrix.columns:
			continue
		total_length += 1

		if aa1 != aa2:
			aa_distance += 1

	return aa_distance/total_length

#calculate protein distance based on PAM matrices and phylogenetic trees
def phylo_distance(seq1: Seq, seq2: Seq, method = "dayhoff"):
	#remove alignment gaps to be consistent with other distance methods in this script
	trimmed_seq1 = MutableSeq("")
	trimmed_seq2 = MutableSeq("")
	for index,aa1 in enumerate(str(seq1)):
		aa2 = seq2[index]
		if(aa1 == "-" or aa2 == "-"):
			continue
		trimmed_seq1.append(aa1)
		trimmed_seq2.append(aa2)

	#Create distance matrix
	seqrec1 = SeqRecord(trimmed_seq1, id="a1", name="a1", )
	seqrec2 = SeqRecord(trimmed_seq2, id="a2", name="a2")
	input_aln = AlignIO.MultipleSeqAlignment([seqrec1, seqrec2])
	calculator = DistanceCalculator(method)
	dist_matrix = calculator.get_distance(input_aln)

	return dist_matrix['a1','a2']

def get_hla_locus(name: str):
	parts = name.replace("HLA-","").split("*")
	return parts[0]


def data_path(name: str):
	this_dir, this_filename = os.path.split(__file__)
	data_dir = os.path.join(this_dir, "data")
	return os.path.join(data_dir, name)


def main(argv):
	parser = OptionParser(usage="usage: %prog [options] HLA-ALLELE1 HLA-ALLELE2", description="Calculates HLA diversity metrics. Pass allele names according HLA nomenclature.")
	parser.add_option("--whole_protein", action="store_true", default=False, dest="whole_protein", help="Use the whole protein sequence of MHC complex. Otherwise only binding grooves will be considered (Exon 2/3 for MHC class I and exon 2 for class II, respectively).")
	parser.add_option("--batch", dest="filename", help="TSV file with one pair of HLA alleles per row. Use four HLA nomenclature.")

	(options, args) = parser.parse_args()

	if(options.filename != None and not exists(options.filename)):
		raise Exception("File %s does not exist!" % options.filename)
	elif(len(args) > 0 and options.filename != None):
		raise Exception("Use either a batch file or pass two HLA alleles via command line!")
	elif(len(args) != 2 and options.filename == None):
		raise Exception("Specify two HLA alleles!")

	print("#ALLELE1", "ALLELE2", "GRANTHAM", "PDISTANCE", "SANDBERG", "DAYHOFF", "JTT", sep='\t')

	if(len(args) == 2): #input command line alleles
		allele1 = args[0]
		allele2 = args[1]

		locus = get_hla_locus(allele1)
		locus2 = get_hla_locus(allele2)
		if locus != locus2:
			raise Exception("The loci HLA-%s and HLA-%s must be the same" % (locus, locus2))

		#load IMGT/HLA protein aligments
		protein_alignments = AlignIO.read(data_path("%s_prot.msf" % locus), "msf")

		seq1 = get_protein_sequence(allele1, protein_alignments, options.whole_protein)
		seq2 = get_protein_sequence(allele2, protein_alignments, options.whole_protein)

		print(allele1, allele2, round(grantham_distance(seq1, seq2), 5), round(p_distance(seq1,seq2), 5), round(sandberg_distance(seq1,seq2), 5), round(phylo_distance(seq1,seq2, "dayhoff"), 5), round(phylo_distance(seq1,seq2, "jones"), 5),  sep='\t')
	else: #input batch alleles
		file = open(options.filename, "r")

		protein_alignments = {
			"A": AlignIO.read(data_path("A_prot.msf"), "msf"),
			"B": AlignIO.read(data_path("B_prot.msf"), "msf"),
			"C": AlignIO.read(data_path("C_prot.msf"), "msf"),
			"DQA1": AlignIO.read(data_path("DQA1_prot.msf"), "msf"),
			"DQB1": AlignIO.read(data_path("DQB1_prot.msf"), "msf"),
			"DRB1": AlignIO.read(data_path("DRB1_prot.msf"), "msf")
		}

		for line in file:
			(allele1, allele2) = line.split('\t')
			allele1 = parse_allele_name(allele1)
			allele2 = parse_allele_name(allele2)

			#Skip entries of different genes
			locus = get_hla_locus(allele1)
			locus2 = get_hla_locus(allele2)
			if locus != locus2:
				print(allele1, allele2, "", "", "", sep="\t")
				continue
			
			seq1 = get_protein_sequence(allele1, protein_alignments[locus], options.whole_protein)
			seq2 = get_protein_sequence(allele2, protein_alignments[locus], options.whole_protein)

			print(allele1, allele2, round(grantham_distance(seq1, seq2), 5), round(p_distance(seq1,seq2), 5), round(sandberg_distance(seq1,seq2), 5), round(phylo_distance(seq1,seq2, "dayhoff"), 5), round(phylo_distance(seq1,seq2, "jones"), 5),  sep='\t')

if __name__ == "__main__":
	main(sys.argv)
