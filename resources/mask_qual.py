from __future__ import division

import sys
import argparse
import os
import operator
import subprocess
import re
from Bio import SeqIO
from Bio.Seq import Seq

def mask_quality(input_fastq,output_fastq,quality):
	# quality filter score, if lower, mask with N
	quality_filter = int(quality)

	# iterate through
	for record in SeqIO.parse(input_fastq,'fastq'):
		# get the scores
		scores = record.letter_annotations["phred_quality"]
		sequence = str(record.seq)
		new_sequence = []
		index = 0
		while index < len(sequence):
			if int(scores[index]) < quality_filter:
				new_sequence.append('N')
			else:
				new_sequence.append(sequence[index])
			index += 1
		# reset sequence
		record.seq = Seq(''.join(new_sequence))

		# write to output file
		SeqIO.write(record, output_fastq, "fastq")


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_fastq',type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--output_fastq',type=argparse.FileType('w'),required=True)
	parser.add_argument('-q','--quality',required=True)
	opts = parser.parse_args(argv)
	mask_quality(opts.input_fastq,opts.output_fastq,opts.quality)
 
if __name__ == '__main__':
	main(sys.argv[1:])