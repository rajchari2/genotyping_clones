# script to calculation mutation rate from SAM file
# Use the FLAG tag, MD tag and CIGAR string for calling
# FLAG tag: Column 2 (make sure this is 0)
# Map start: Column 4 (Should be position 1)
# CIGAR: Column 6
# MD tag: Column 13

from __future__ import division

import sys
import argparse
import os
import operator
import subprocess
import re
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

from collections import defaultdict
from collections import Counter
from operator import itemgetter

def genotype_sample(bam_file,variant_list,output_file):
	# read the sorted bam file
	samfile = pysam.AlignmentFile(bam_file,"rb")
	sample_name = bam_file.name

	# parse variant list
	var_list = variant_list.split(',')
	proportion_list = []
	for var in var_list:
		# reference base
		ref = var[0]
		variant = var[-1]
		position = int(var[1:-1])

		# coverage goes A, C, G, T
		coverage_data = samfile.count_coverage(start=position-1, end=position)

		# count for each
		a_count = int(coverage_data[0][0])
		c_count = int(coverage_data[1][0])
		g_count = int(coverage_data[2][0])
		t_count = int(coverage_data[3][0])
		total = a_count + c_count + g_count + t_count

		# check variant proportion
		value = 0
		if variant=='A':
			value = a_count
		elif variant=='C':
			value = c_count
		elif variant=='G':
			value = g_count
		else:
			value = t_count
		proportion = value / total
		proportion_list.append(proportion)
	
	# calculate the average
	average = sum(proportion_list) / len(proportion_list)
	output_file.write(sample_name  + '\t' + str(average) + '\n')

	# close handle
	bam_file.close()
	


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--sample_bam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-v','--variant_list',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	genotype_sample(opts.sample_bam_file,opts.variant_list,opts.output_file)

if __name__ == '__main__':
	main(sys.argv[1:])
