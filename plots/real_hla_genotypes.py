import glob
import os
import sys


def get_real_pairs():
	p_dir = "/mnt/storage2/GRCh38/projects/diagnostic/SomaticAndTreatment/"
	som_dirs = glob.glob(p_dir + "Somatic_*-*")

	o_handle = open("real_allele_pairs.txt", "w")

	for dir in som_dirs:
		(tid, nid) = os.path.basename(dir).replace("Somatic_", "").split("-")
		n_hla_file = os.path.join(p_dir, "Sample_" + nid, nid + "_hla_genotyper.tsv")
		if not os.path.exists(n_hla_file):
			continue
		handle = open(n_hla_file, "r")
		for line in handle:
			parts = line.strip().split("\t")
			(gene,a1,a2, pval, qual) = parts[3:8]
			if not gene in ["HLA-A", "HLA-B", "HLA-C"]:
				continue
			
			if qual != "Pass":
				continue

			o_handle.write("%s\t%s\n" % (a1,a2))
		handle.close()
	
	o_handle.close()


def main(argv):
	get_real_pairs()


if __name__ == "__main__":
	main(sys.argv)