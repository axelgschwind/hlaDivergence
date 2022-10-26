all: download_alignments download_common_alleles


#Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE: IPD-IMGT/HLA Database. Nucleic Acids Research (2020), 48:D948-55
#Robinson J, Malik A, Parham P, Bodmer JG, Marsh SGE: IMGT/HLA - a sequence database for the human major histocompatibility complex Tissue Antigens (2000), 55:280-287
download_alignments:
	wget ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/A_prot.msf
	wget ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/B_prot.msf
	wget ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/C_prot.msf
	wget ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/DRB1_prot.msf
	wget ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/DQA1_prot.msf
	wget ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/DQB1_prot.msf

#Hurley CK et al: Common, intermediate and well-documented HLA alleles in world populations: CIWD version 3.0.0. HLA. 2020 Jun;95(6):516-531. doi: 10.1111/tan.13811. Epub 2020 Jan 31. PMID: 31970929; PMCID: PMC7317522.
download_common_alleles:
	wget -O ciwd3.xlsx https://s3.eu-central-1.amazonaws.com/ihiw.website.data/CIWD-3.0/CIWD_Palleles-all_loci-2020320-ihws-website.xlsx
	echo "\33[41m You have to parse the excel file manually into a txt file which contains one common allele name per row.\n \
After downloading, the excel sheet has to be filtered to common alleles (Column 'Highest Frequency' = 'C') and the allele names have to be copied to the file common_alleles.txt manually.\n \
Every row must contain a common allele name\n\033[m"
	python ../group_common_alleles.py > grouped_common_alleles.tsv