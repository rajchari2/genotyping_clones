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

def process_cigar(cigar_tuple,md_tag,aligned_pairs,query_qualities):
	alterations = []
	pos = 1
	# if the cigar tuple is only 1 entry check if it has a SNP
	if len(cigar_tuple)==1:
		for seq_pos in aligned_pairs:
			pos = int(seq_pos[0]) + 1
			base = seq_pos[2]
			if base.islower() and int(query_qualities[pos-1]) >= 20:
				alt = str(pos) + ':' + str(pos) + ':0:SNP'
				alterations.append(alt)
	else:
		for entry in cigar_tuple:
			# match
			length = int(entry[1])
			if int(entry[0])==0:
				pos = pos + int(entry[1])
			elif int(entry[0])==1:
				alt = str(pos) + ':' + str(pos) + ':' + str(length) + ':Ins'
				alterations.append(alt)
			elif int(entry[0])==2:
				end_pos = pos + int(entry[1])
				alt = str(pos) + ':' + str(end_pos) + ':' + str(length) + ':Del'
				pos = end_pos
				alterations.append(alt)
	return alterations

def find_target_indices(target_site,reference_file,coordinates):
	# store reference sequence
	# store reference sequence
	refDB = defaultdict(str)
	targetDB = defaultdict(str)
	gene = ''
	if coordinates != 'full':
		chrom,pos = coordinates.split(':')
		start,end = pos.split('-')
	for record in SeqIO.parse(reference_file,'fasta'):
		if coordinates=='full':
			refDB[str(record.id)] = str(record.seq)
			targetDB[str(record.id)] = target_site
			gene = str(record.id)
			break
		elif chrom==str(record.id):
			refDB[str(record.id)] = str(record.seq)[int(start):int(end)]
			targetDB[str(record.id)] = target_site
			gene = str(record.id)
			break

	# identify the target site start and end
	targetSeq = targetDB[gene].upper()
	refSeq = refDB[gene].upper()

	# identify where in the sequence it is
	tSeq = Seq(targetSeq)
	targetSeqRC = str(tSeq.reverse_complement())

	# find the sequence
	index = refSeq.find(targetSeq)

	if index==-1:
		rcIndex = refSeq.find(targetSeqRC)
		startIndex = rcIndex
		endIndex = startIndex + len(targetSeqRC)
	else:
		startIndex = index	
		endIndex = startIndex + len(targetSeq)

	return startIndex,endIndex,gene



def genotype_sample(bam_file,control_bam_file,expt_type,variant_list,target_site,reference_file,coordinates,output_file):
	# read the sorted bam file
	samfile = pysam.AlignmentFile(bam_file,"rb")
	# go through control file to determine "false" mutations
	# control bam file
	ctrl_bam_file = 'processed/' + control_bam_file + '_bwamem_sorted.bam'
	control_sam_file = pysam.AlignmentFile(ctrl_bam_file,'rb')
	sample_name = bam_file.name
	sample_name = sample_name.replace('processed/','')
	sample_name = sample_name.replace('_bwamem_sorted.bam','')

	if expt_type=='HDR':

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
			coverage_data_ctrl = control_sam_file.count_coverage(start=position-1, end=position)

			# count for each
			a_count = int(coverage_data[0][0])
			c_count = int(coverage_data[1][0])
			g_count = int(coverage_data[2][0])
			t_count = int(coverage_data[3][0])
			total = a_count + c_count + g_count + t_count

			# count for control
			total = int(coverage_data_ctrl[0][0]) + int(coverage_data_ctrl[1][0]) + int(coverage_data_ctrl[2][0]) + int(coverage_data_ctrl[3][0])

			# check variant proportion
			value = 0
			if variant=='A':
				value = a_count
				ctrl_value = int(coverage_data_ctrl[0][0])
			elif variant=='C':
				value = c_count
				ctrl_value = int(coverage_data_ctrl[1][0])
			elif variant=='G':
				value = g_count
				ctrl_value = int(coverage_data_ctrl[2][0])
			else:
				value = t_count
				ctrl_value = int(coverage_data_ctrl[3][0])
			proportion = value / total
			proportion_list.append(proportion)

			ctrl_proportion = ctrl_value / total
			ctrl_proportion_list.append(ctrl_proportion)
		
		# calculate the average
		average = sum(proportion_list) / len(proportion_list)
		ctrl_average = sum(ctrl_proportion_list) / len(ctrl_proportion_list)
		corrected_average = average - ctrl_average
		output_file.write(sample_name  + '\t' + str(corrected_average) + '\n')

	# in case of deletions, just count indels
	else:
		# variables to store information
		sample_read_total = 0
		sample_mutation_count = 0
		control_cigar_strings = []

		# get the start and end indexes of the target in the reference
		target_start,target_end,gene = find_target_indices(target_site,reference_file,coordinates)

		# get reference position
		if coordinates=='full':
			ref_pos = 0
		else:
			chrom,pos = coordinates.split(':')
			start,end = pos.split('-')
			ref_pos = int(start)-1

		# write header in output file
		output_file.write('Sample\tNHEJ_Mutation_Rate\n')
		for read in control_sam_file.fetch():
			if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==ref_pos:
				# check if the read has a mutation that involves the target sequence
				md_tag = read.get_tag('MD')
				ctrl_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities)
				unique_entry = md_tag + '_' + read.cigarstring
				# if alteration is non zero
				valid_alteration = False
				if len(ctrl_alterations) > 0:
					for alt in ctrl_alterations:
						[start,end,size,alt_type] = alt.split(':')
						if ((int(start) >= target_start and int(start) <= target_end) or (int(end) >= target_start and int(end) <= target_end)):
							valid_alteration = True
				# check if valid alteration
				if valid_alteration==True and unique_entry not in control_cigar_strings:
					# add to the control
					control_cigar_strings.append(unique_entry)

		# go through sample BAM file
		for read in samfile.fetch():
			if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==ref_pos:
				sample_read_total += 1
				md_tag = read.get_tag('MD')
				sample_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities)
				unique_entry = md_tag + '_' + read.cigarstring
				# only count if the tuple is not in the control
				if len(sample_alterations) > 0 and unique_entry not in control_cigar_strings:
					# if alteration is non zero
					valid_alteration = False
					for alt in sample_alterations:
						[start,end,size,alt_type] = alt.split(':')
						if ((int(start) >= target_start and int(start) <= target_end) or (int(end) >= target_start and int(end) <= target_end)) and int(size) % 3 != 0:
							valid_alteration = True
					if valid_alteration==True:
						sample_mutation_count += 1
		if sample_read_total >= 100:
			rate = (sample_mutation_count / sample_read_total) * 100
		else:
			rate = 'N/A'

		# write to output
		output_file.write(sample_name + '\t' + str(rate) + '\n')

	# close handle
	bam_file.close()
	output_file.close()


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--sample_bam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-c','--control_bam_file',required=True)
	parser.add_argument('-e','--experiment_type',required=True)
	parser.add_argument('-v','--variant_list',required=True)
	parser.add_argument('-t','--target_site',required=True)
	parser.add_argument('-r','--reference_file',required=True)
	parser.add_argument('-p','--coordinates',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	genotype_sample(opts.sample_bam_file, opts.control_bam_file, opts.experiment_type, opts.variant_list, opts.target_site, opts.reference_file, opts.coordinates, opts.output_file)

if __name__ == '__main__':
	main(sys.argv[1:])
