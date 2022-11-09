import matplotlib.pyplot as plt
import pandas as pd
import itertools
import sys
from optparse import OptionParser
from scipy.stats import norm, pearsonr, spearmanr
import numpy as np


data_ = pd.DataFrame()


def scatterplot_pairwise_metrics(loci = ["all", "A", "B", "C", "DQB1", "DRB1"]):
	global data_
	metrics = [metric for metric in itertools.combinations(["PDISTANCE", "JTT", "DAYHOFF", "GRANTHAM", "SANDBERG"], 2) ]

	for locus in loci:
		#filter data per HLA locus
		data_per_locus = data_.loc[data_["ALLELE1"].str.contains("%s\*" % locus)]
		if locus == "all":
			data_per_locus = data_

		fig = plt.figure(figsize=(8,3.7)) # determines aspect ratio of the whole plot

		fig.suptitle("HLA-%s allele pairwise divergence scatter plots" % locus)


		for i, (metric1, metric2) in enumerate(metrics, start=1):
			ax1 = fig.add_subplot(2,5,i)

			data_per_locus.plot( x=metric1, y=metric2, kind="scatter", s=1, ax=ax1)
			
			#set 1:1 aspect ratio
			ax1.set_aspect( 1./ax1.get_data_ratio())

			ax1.set_xlabel(metric1, fontsize="xx-small")
			ax1.set_ylabel(metric2, fontsize="xx-small")
			ax1.tick_params(axis='both', which='major', labelsize="xx-small")

			#Calculate Pearson's r correlation coefficient
			(pearson_r, pearson_p) = pearsonr(data_per_locus[metric1], data_per_locus[metric2])

			ax1.text(0.1, 0.95, "$r=%.2f$" % pearson_r, transform=ax1.transAxes, fontsize="xx-small", verticalalignment='top')

		fig.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
		fig.savefig("HLA-%s_pairwise_correlation.png" % locus, dpi=300)

		plt.close()

def histogram(data: pd.DataFrame, prefix = "", loci = ["all", "A", "B", "C", "DQB1", "DRB1"], metrics = ["PDISTANCE", "JTT", "DAYHOFF", "GRANTHAM", "SANDBERG"]):

	for locus in loci:
		data_per_locus = data.loc[data["ALLELE1"].str.contains("%s\*" % locus)]
		if locus == "all":
			data_per_locus = data

		fig = plt.figure()

		fig.subplots_adjust(hspace=0.5, wspace=0.5)
		fig.suptitle("HLA-%s allele divergence" % locus)

		for i, metric in enumerate(metrics, start=1):
			ax = fig.add_subplot(2,3,i)

			ax.set_xlabel(metric)
			ax.set_ylabel(r'$p(%s)$' % metric)

			ax.hist(data_per_locus[metric], bins = 30, density=True, edgecolor='black', linewidth=1)

			mean, sd = norm.fit(data_per_locus[metric])
			xmin,xmax = plt.xlim()
			x = np.linspace(xmin, xmax, 100)

			p = norm.pdf(x, mean, sd)
			
			# place a text box in upper left in axes coords
			ax.text(0.6, 0.95, "$\mu=%.2f$" % mean, transform=ax.transAxes, fontsize="small", verticalalignment='top')

			ax.plot(x, p, color="red", linewidth=1)

		fig.savefig("%s%s_divergence_histograms.png" % (prefix,locus), dpi=300)


def main(argv):
	global data_
	parser = OptionParser(usage="usage: %prog [options] HLA locus", description="Several plots for HLA diversity scores.")
	parser.add_option("--distances", dest="filename", help="TSV file with HLA distances as produces by main.py.", default="common_distances_binding_grooves.tsv")
	parser.add_option("--prefix", dest="prefix", help="Prefix for output plots", default= "")
	(options, args) = parser.parse_args()

	print(options.prefix)

	data_ = pd.read_csv(options.filename, sep="\t", header=0)

	data_.rename(columns={"#ALLELE1": "ALLELE1"}, inplace=True)

	if len(args) > 0:
		scatterplot_pairwise_metrics(args)
		histogram(data_, options.prefix, args)
	else:
		scatterplot_pairwise_metrics()
		histogram(data_, options.prefix)

if __name__ == "__main__":
	main(sys.argv)