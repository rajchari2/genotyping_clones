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

def process_cigar(cigar_tuple,md_tag,aligned_pairs,query_qualities,query_alignment_sequence,query_sequence):
	alterations = []
	pos = 1
	# if the cigar tuple is only 1 entry check if it has a SNP
	if len(cigar_tuple)==1:
		index = 0
		while index < len(aligned_pairs):
			current_tuple = aligned_pairs[index]
			pos = int(current_tuple[0]) + 1
			base = current_tuple[2]
			next_tuple = ()
			start_index = index

			snp_found = False
			while base.islower() and int(query_qualities[pos-1]) >= 20 and index < len(aligned_pairs)-1:
				snp_found = True
				index = index + 1
				next_tuple = aligned_pairs[index]
				pos = int(next_tuple[0]) + 1
				base = next_tuple[2]
			
			if snp_found==True:
				# reset the index by 1
				index = index - 1
				ref_seq = ''
				start_pos = int(aligned_pairs[start_index][0]) + 1 
				end_pos = int(aligned_pairs[index][0]) + 1
				start_chrom_pos = aligned_pairs[start_index][1]+1
				final_seq = query_alignment_sequence[aligned_pairs[start_index][0]:aligned_pairs[index][0]+1]
				end_chrom_pos = aligned_pairs[index][1] + 1
				size = end_pos - start_pos + 1
				while start_index <= index:
					my_tuple = aligned_pairs[start_index]
					ref_seq = ref_seq + my_tuple[2].upper()
					start_index += 1
				alt = str(start_pos) + ':' + str(start_pos) + ':0:' + str(start_chrom_pos) + '-' + str(end_chrom_pos) + ':SNV-' + ref_seq + '>' + final_seq
				alterations.append(alt)
				#print(alt)
			index += 1
	else:
		for entry in cigar_tuple:
			if int(entry[0])==0:
				pos = pos + int(entry[1])
			elif int(entry[0])==1:
				index = pos-1
				ins_start = aligned_pairs[index-1][1]
				ins_seq = query_alignment_sequence[index-1:index-1+int(entry[1])]
				alt = str(pos) + ':' + str(pos) + ':' + str(entry[1]) + ':' + str(ins_start+1) + '-' + str(ins_start+1) + ':Ins-' + ins_seq
				alterations.append(alt)
				pos = pos + int(entry[1])
			elif int(entry[0])==2:
				index = pos - 1
				current_tuple = aligned_pairs[index]
				del_start = current_tuple[1]
				del_end = current_tuple[1]
				del_seq = ''
				while index < len(aligned_pairs) and current_tuple[0]==None:
					del_seq = del_seq + current_tuple[2]
					del_end = current_tuple[1]
					index += 1
					current_tuple = aligned_pairs[index]					
				end_pos = pos + int(entry[1])
				# print(end_pos)
				alt =  str(pos) + ':' + str(end_pos) + ':' + str(entry[1]) + ':' + str(del_start+1) + '-' + str(del_end+1) + ':Del-' + del_seq
				pos = end_pos
				alterations.append(alt)
	return alterations

def find_target_indices(target_site,reference_file,coordinates):
	# store reference sequence
	refDB = defaultdict(str)
	targetDB = defaultdict(str)
	gene = ''
	chrom,pos = coordinates.split(':')
	start,end = pos.split('-')
	for record in SeqIO.parse(reference_file,'fasta'):
		if chrom==str(record.id):
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


def genotype_sample(sample_sam_file,control,expt_type,variant_list,target_site,reference_file,coordinates,output_file):

	# variables to store information
	sample_read_total = 0
	control_cigar_strings = []

	# control bam file
	control_sam_file = 'processed/' + control + '_bwamem_sorted.bam'

	# open file handles
	ctrl_sam = pysam.AlignmentFile(control_sam_file,'rb')
	samfile = pysam.AlignmentFile(sample_sam_file,"rb")

	# get the start and end indexes of the target in the reference
	target_start,target_end,gene = find_target_indices(target_site,reference_file,coordinates)

	# get reference position
	chrom,pos = coordinates.split(':')
	start,end = pos.split('-')
	ref_pos = int(start)-1

	# write header in output file
	output_file.write('Sample\tTarget_Site\tTotal_Read_Count\tNHEJ_Mutation_Rate\tOOF_Mutation_Rate\tHDR_Rate\tAlterations\tCheckForLargeIndel\n')

	# go through control file to determine "false" mutations		
	ctrl_total_count = 0
	ctrl_mut_freq = defaultdict(int)
	for read in ctrl_sam.fetch():
		if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==ref_pos:
			# check if the read has a mutation that involves the target sequence
			md_tag = read.get_tag('MD')
			ctrl_total_count += 1
			ctrl_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities,read.query_alignment_sequence, read.query_sequence)
			unique_entry = md_tag + '_' + read.cigarstring
			# if alteration is non zero
			valid_alteration = False
			if len(ctrl_alterations) > 0:
				for alt in ctrl_alterations:
					[start,end,indel_len,indel_coords,alt_type] = alt.split(':')
					if ((int(start) >= target_start-3 and int(start) <= target_end+3) or (int(end) >= target_start-3 and int(end) <= target_end+3)):
						valid_alteration = True
			# check if valid alteration
			if valid_alteration==True:
				if unique_entry not in ctrl_mut_freq:
					ctrl_mut_freq[unique_entry] = 1
					control_cigar_strings.append(unique_entry)
				else:
					ctrl_mut_freq[unique_entry] += 1

	# go through actual sample file
	unique_alts = defaultdict(list)
	cigar_counts = defaultdict(int)
	is_oof = defaultdict(str)
	for read in samfile.fetch():
		if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==ref_pos:
			sample_read_total += 1
			md_tag = read.get_tag('MD')
			sample_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities,read.query_alignment_sequence, read.query_sequence)
			unique_entry = md_tag + '_' + read.cigarstring

			# count by cigar string
			if read.cigarstring not in cigar_counts:
				cigar_counts[read.cigarstring] = 1
			else:
				cigar_counts[read.cigarstring] += 1

			# only count if the tuple is not in the control at a high frequency
			entry_freq = 0
			if unique_entry in ctrl_mut_freq and ctrl_total_count > 0:
				entry_freq = ctrl_mut_freq[unique_entry] / ctrl_total_count

			# if it's less than 10% in the control, we can still count it
			if len(sample_alterations) > 0 and entry_freq < 0.1:
				# if alteration is non zero
				valid_alteration = False
				for alt in sample_alterations:
					[start,end,indel_len,indel_coords,alt_type] = alt.split(':')
					if ((int(start) >= target_start-3 and int(start) <= target_end+3) or (int(end) >= target_start-3 and int(end) <= target_end+3)):
						valid_alteration = True
						final_alt = chrom + ':' + indel_coords + '&' + alt_type
						unique_alts[read.cigarstring].append(final_alt)
						if int(indel_len) % 3 != 0:
							is_oof[final_alt] = 'Yes'
						break

	# list of alterations
	alt_line = []

	# calculate mutation and OOF and correct for spurious alterations
	if sample_read_total < 20:
		nhej_rate = 'N/A'
		oof_rate = 'N/A'
		sample_read_total_corrected = sample_read_total
	else:
		# make copies of the variable so that they can be adjusted
		sample_read_total_corrected = 0
		sample_mutation_count_corrected = 0
		sample_oof_count_corrected = 0

		# go through the most prevalent cigars
		valid_cigars = []
		for cigar in cigar_counts:
			if cigar_counts[cigar] / sample_read_total >= 0.2:
				sample_read_total_corrected += cigar_counts[cigar]
				valid_cigars.append(cigar)

		# go through the prevalent cigars and determine alteration
		alt_freq = defaultdict(int)
		for cigar in valid_cigars:
			for alt in unique_alts[cigar]:
				sample_mutation_count_corrected += 1
				if alt not in alt_freq:
					alt_freq[alt] = 1
				else:
					alt_freq[alt] += 1
				if alt in is_oof:
					sample_oof_count_corrected += 1

		# calculate rates based on corrected numbers
		if sample_read_total_corrected > 0:
			nhej_rate = (sample_mutation_count_corrected / sample_read_total_corrected) * 100
			oof_rate = (sample_oof_count_corrected / sample_read_total_corrected) * 100

			# get the list of alterations
			sorted_alts = sorted(alt_freq.items(), key=lambda kv: kv[1], reverse=True)
			for u_alt in sorted_alts:
				proportion = round((alt_freq[u_alt[0]] / sample_read_total_corrected * 100),2)
				if proportion >= 20:
					alt = u_alt[0] + '&' + str(proportion)
					alt_line.append(alt)
		else:
			nhej_rate = 'N/A'
			oof_rate = 'N/A'

	hdr_rate = 'N/A'
	if expt_type=='HDR':

		# parse variant list
		proportion_list = []
		ctrl_proportion_list = []
		for var in variant_list:
			# split them into their parts
			alt,mut_base,chrom,position = var.split(',')

			# first deal with the easiest, SNP
			output_file.write('Sample\tProportions\n') 
			if alt=='SNV':
				# coverage goes A, C, G, T
				coverage_data = samfile.count_coverage(chrom,start=int(position)-1, end=int(position))
				coverage_data_ctrl = ctrl_sam.count_coverage(chrom,start=int(position)-1, end=int(position))

				# count for each
				a_count = int(coverage_data[0][0])
				c_count = int(coverage_data[1][0])
				g_count = int(coverage_data[2][0])
				t_count = int(coverage_data[3][0])
				total_sample = a_count + c_count + g_count + t_count

				# count for control
				total_control = int(coverage_data_ctrl[0][0]) + int(coverage_data_ctrl[1][0]) + int(coverage_data_ctrl[2][0]) + int(coverage_data_ctrl[3][0])

				# check variant proportion
				value = 0
				if mut_base=='A':
					value = a_count
					ctrl_value = int(coverage_data_ctrl[0][0])
				elif mut_base=='C':
					value = c_count
					ctrl_value = int(coverage_data_ctrl[1][0])
				elif mut_base=='G':
					value = g_count
					ctrl_value = int(coverage_data_ctrl[2][0])
				else:
					value = t_count
					ctrl_value = int(coverage_data_ctrl[3][0])
				proportion = value / total_sample
				proportion_list.append(proportion)

				ctrl_proportion = ctrl_value / total_control
				ctrl_proportion_list.append(ctrl_proportion)
		
		# calculate the average
		hdr_rate = str(sum(proportion_list) / len(proportion_list) * 100)
		average_control = sum(ctrl_proportion_list) / len(ctrl_proportion_list) * 100

	# write to output file
	large_indel = 'No'
	if len(alt_line)==1:
		large_indel = 'Yes'

	# update the sample name
	sample_name = sample_sam_file.name.replace('processed/','')
	sample_name = sample_name.replace('_bwamem_sorted.bam','')

	# write to output file
	output_file.write(sample_name + '\t' + target_site + '\t' + str(sample_read_total_corrected) + '\t' + str(nhej_rate) + '\t' + str(oof_rate) + '\t' + hdr_rate + '\t' + ';'.join(alt_line) + '\t' + large_indel + '\n')

	# close handle
	samfile.close()
	output_file.close()


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--sample_bam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-c','--control_bam_file',required=True)
	parser.add_argument('-e','--experiment_type',required=True)
	parser.add_argument('-v','--variant_list',nargs='+',required=True)
	parser.add_argument('-t','--target_site',required=True)
	parser.add_argument('-r','--reference_file',required=True)
	parser.add_argument('-p','--coordinates',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	genotype_sample(opts.sample_bam_file, opts.control_bam_file, opts.experiment_type, opts.variant_list, opts.target_site, opts.reference_file, opts.coordinates, opts.output_file)

if __name__ == '__main__':
	main(sys.argv[1:])
