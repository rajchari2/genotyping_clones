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

	# go through mutation file list
	for infile in input_file_list:
		ifile = open(infile,'r')
		for line in ifile:
			line = line.rstrip('\r\n')
			toml_file.write(line + '\n')
		ifile.close()

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