# script store all alterations not observed in WT or Plasmid library samples
# updated to not count spurious SNPs
# updated to not count spurious indels too
# updated to not count SNPs at all (only look at insertions and deletions, even if they happen once)
# edited: March 2018

from __future__ import division

import sys
import argparse
import os
import operator
import subprocess
import re
from Bio import SeqIO
from Bio.Seq import Seq

from collections import defaultdict
from collections import Counter
from operator import itemgetter

def aggregate_files(input_file_list,toml_file):
	# data structure
	plate_data = defaultdict(dict)

	# go through mutation file list
	for infile in input_file_list:
		ifile = open(infile,'r')
		line_count = 0
		for line in ifile:
			if line_count > 0:
				line = line.rstrip('\r\n')
				parts = line.split('\t')
				# parse sample name
				sample_name_parts = parts[0].split('-')
				plate = sample_name_parts[0] + sample_name_parts[2]
				well = sample_name_parts[1]
				# initialize
				if plate not in plate_data:
					plate_data[plate] = defaultdict(float)
				if parts[1]=='N/A':
					value = 0.0
				else:
					value = float(parts[1])
				plate_data[plate][well] = value
			line_count += 1
		ifile.close()

	# write to the toml file
	for plate in plate_data:
		toml_file.write('[plate.' +  plate + ']' + '\n')
		for well in plate_data[plate]:
			toml_file.write('  well.' + well + '.genotype=' + str(plate_data[plate][well]) + '\n')

	# close file handles
	toml_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file_list',nargs='+',required=True)
	parser.add_argument('-o','--toml_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	aggregate_files(opts.input_file_list,opts.toml_file)

if __name__ == '__main__':
	main(sys.argv[1:])