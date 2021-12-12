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

def zygosity(value):
	zygo = ''
	if value == 'N/A':
		zygo = 'NC'
	elif value >= 20 and value <= 80:
		zygo = 'Heterozygous'
	elif value <= 5:
		zygo = 'WT'
	elif value >= 90:
		zygo = 'Homozygous'
	return zygo

def aggregate_files(input_file_list,expt_type,toml_file,aggregated_file,clone_list):
	# data structure
	plate_data = defaultdict(dict)
	clone_list.write('Sample\tTarget_Site\tTotal_Read_Count\tNHEJ_Mutation_Rate\tOOF_Mutation_Rate\tHDR_Rate\tAlterations\tCheckForLargeIndel\tGenotype\n')
	aggregated_file.write('Sample\tTarget_Site\tTotal_Read_Count\tNHEJ_Mutation_Rate\tOOF_Mutation_Rate\tHDR_Rate\tAlterations\tCheckForLargeIndel\tGenotype\n')

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
				plate = sample_name_parts[0] + '-' + sample_name_parts[2]
				well = sample_name_parts[1]

				# get values
				nhej_rate = parts[3]
				oof_rate = parts[4]
				hdr_rate = parts[5]

				# initialize
				if plate not in plate_data:
					plate_data[plate] = defaultdict(str)

				# possible genotypes for KO: WT, heterozygous, homozygous-NHEJ, homozygous-OOF, NC
				# possible genotypes for KI: WT/WT, WT/KO, KO/KO, KI/KI, KI/KO, WT/KI, NC
				value = ''
				if expt_type=='KO':
					if nhej_rate=='N/A' or oof_rate=='N/A':
						value = 'NC'
					elif float(nhej_rate) <= 10:
						value = 'WT'
					elif float(nhej_rate) >= 90 and float(oof_rate) < 90:
						value = 'Homozygous-NHEJ'
					elif float(nhej_rate) >= 90 and float(oof_rate) >= 90:
						value = 'Homozygous-OOF'
					else:
						value = 'Heterozygous'
				else:
					if nhej_rate=='N/A' or oof_rate=='N/A' or hdr_rate=='N/A':
						value = 'NC'
					else:
						zygo_nhej = zygosity(nhej_rate)
						zygo_hdr =  zygosity(hdr_rate)
						zygo_oof = zygosity(oof_rate)
						if zygo_nhej=='WT' and zygo_hdr=='WT':
							value = 'WT/WT'
						elif zygo_nhej=='Heterozygous' and zygo_oof=='WT' and zygo_hdr=='Heterozygous':
							value = 'WT/KI'
						elif zygo_nhej=='Homozygous' and zygo_hdr=='Homozygous':
							value = 'KI/KI'
						elif zygo_nhej=='Homozygous' and zygo_hdr=='Heterozygous':
							value = 'KI/KO'
						elif zygo_nhej=='Homozygous' and zygo_hdr=='WT':
							value = 'KO/KO'
						elif zygo_nhej=='Heterozygous' and zygo_hdr=='WT':
							value = 'WT/KO'

				# record the value
				plate_data[plate][well] = value

				# write out clone list
				if (expt_type=='KO' and (value=='Homozygous-NHEJ' or value=='Homozygous-OOF')) or (expt_type=='HDR' and 'KI' in value):
					clone_list.write(line + '\t' + value + '\n')

				# write to the aggregated output
				aggregated_file.write(line + '\t' + value +'\n')
			line_count += 1
		ifile.close()
	clone_list.close()

	# write to the toml file
	for plate in plate_data:
		toml_file.write('[plate.' +  plate + ']' + '\n')
		for well in plate_data[plate]:
			toml_file.write('  well.' + well + '.genotype=\'' + plate_data[plate][well] + '\'\n')

	# close file handles
	toml_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file_list',nargs='+',required=True)
	parser.add_argument('-e','--expt_type',required=True)
	parser.add_argument('-c','--clone_list',type=argparse.FileType('w'),required=True)
	parser.add_argument('-a','--aggregated_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-o','--toml_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	aggregate_files(opts.input_file_list,opts.expt_type,opts.toml_file,opts.aggregated_file,opts.clone_list)

if __name__ == '__main__':
	main(sys.argv[1:])