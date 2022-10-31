import matplotlib.pyplot as plt
import pandas as pd
import itertools
import sys


def scatterplot_pairwise_metrics(loci = ["all", "A", "B", "C", "DQA1", "DQB1", "DRB1"]):
	data = pd.read_csv("distances.tsv", sep="\t", header=0)
	data.rename(columns={"#ALLELE1": "ALLELE1"}, inplace=True)
	metrics = [metric for metric in itertools.combinations(["PDISTANCE", "JTT", "DAYHOFF", "GRANTHAM", "SANDBERG",], 2) ]

	for locus in loci:
		data_per_locus = data.loc[data["ALLELE1"].str.contains("%s\*" % locus)]
		if locus == "all":
			data_per_locus = data
		
		for (metric1, metric2) in metrics:
			ax = data_per_locus.plot.scatter( x=metric1, y=metric2 )
			ax.set_xlabel(metric1)
			ax.set_ylabel(metric2)

			plt.savefig("%s_%s_%s.png" %  (locus, metric1, metric2) )
			plt.close()

def main(argv):
	loci = []
	if len(argv) > 1:
		for i in range(1, len(argv)):
			loci.append(argv[i])
		scatterplot_pairwise_metrics(loci)
	else:
		scatterplot_pairwise_metrics()

	print(loci)
if __name__ == "__main__":
	main(sys.argv)