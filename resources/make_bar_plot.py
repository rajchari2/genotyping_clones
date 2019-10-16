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

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

def make_indel_rate_plot(input_file,output_file):
	# read in the data
	x_data_labels = []
	y_data = []
	line_count = 0
	mapped_data = defaultdict(float)
	for line in input_file:
		if line_count > 0:
			line = line.rstrip('\r\n')
			parts = line.split('\t')
			mapped_data[parts[0]] = float(parts[5])
		line_count += 1
	input_file.close()

	# put the data into arrays
	for key in sorted(mapped_data.keys()):
		x_data_labels.append(key)
		y_data.append(mapped_data[key])

	# number of bars
	num_bars = len(x_data_labels)
	x_locs = range(num_bars)

	# add axis labels
	fig, ax = plt.subplots()
	plt.bar(x_locs,y_data,color='b')
	plt.xlabel('Guide RNA Candidate')
	plt.ylabel('NHEJ-induced mutation rate (%)')
	ax.set_xticks(x_locs)
	ax.set_xticklabels(x_data_labels,rotation=90)

	# save to file
	plt.savefig(output_file)

def make_diversity_plot(diversity_summary_file,diversity_graph_output):
	# read in the data
	x_data_labels = []
	y_data = []
	line_count = 0
	deletion_data = defaultdict(float)
	insertion_data = defaultdict(float)
	total_data = defaultdict(float)
	mut_data = defaultdict(str)
	for line in diversity_summary_file:
		if line_count > 0:
			line = line.rstrip('\r\n')
			parts = line.split('\t')
			del_rate = 0.0
			ins_rate = 0.0
			total_rate = 0.0
			mut_data[parts[0]] = parts[3]
			if int(parts[3]) != 0:
				del_rate = int(parts[4]) / int(parts[3])
				ins_rate = int(parts[5]) / int(parts[3])
				total_rate = int(parts[6]) / int(parts[3])
			# set the values
			deletion_data[parts[0]] = del_rate
			insertion_data[parts[0]] = ins_rate
			total_data[parts[0]] = total_rate
		line_count += 1
	diversity_summary_file.close()

	# add into lists sorted by key
	x_labels = []
	del_array = []
	ins_array = []
	total_array = []
	for key in sorted(deletion_data.keys()):
		name = key + '_' + mut_data[key]
		x_labels.append(name)
		del_array.append(deletion_data[key])
		ins_array.append(insertion_data[key])
		total_array.append(total_data[key])

	# plotting
	bar_width = 1

	r1 = np.arange(len(x_labels)) * 4
	r2 = r1 + bar_width
	r3 = r2 + bar_width
	
	# add axis labels
	fig, ax = plt.subplots()

	plt.bar(r1, del_array, color='k', width=bar_width, label='Deletions')
	plt.bar(r2, ins_array, color='m', width=bar_width, label='Insertions')
	plt.bar(r3, total_array, color='b', width=bar_width, label='Total')

	plt.xlabel('Guide RNA Candidate - # of mutated reads')
	plt.ylabel('Diversity proportion')
	ax.set_xticks(r1+1)
	ax.set_xticklabels(x_labels,rotation=90)

	# save to file
	ax.autoscale_view()
	ax.legend()
	plt.savefig(diversity_graph_output)

def make_base_editing_plot(input_file,modality,output_file):
	# only going to use position 5 to 12
	# output_file.write('File\tGene\tTotalCount\tPos1\tPos2\tPos3\tPos4\tPos5\tPos6\tPos7\tPos8\tPos9\tPos10\tPos11\tPos12\tPos13\tPos14\tPos15\tPos16\tPos17\tPos18\tPos19\tPos20\n')
	# data structure to hold data
	data_table = defaultdict(list)
	legend_list = []
	# go through each summary
	line_count = 0
	for line in input_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		if line_count > 0:
			depth = int(parts[2])
			index_to_start = 3
			while index_to_start < len(parts):
				value = int(parts[index_to_start]) / depth
				data_table[parts[0]].append(value)
				index_to_start += 1
		else:
			x_data_labels = parts[3:]
		line_count += 1

	# bar locations
	ind = np.arange(0,3*len(x_data_labels),3)
	width = 0.3
	
	# add axis labels
	fig, ax = plt.subplots(figsize=(20,10))
	ax.set_xticks(ind)
	ax.set_xticklabels(x_data_labels,rotation=90)

	# iterate through data
	iteration = 1
	for key in sorted(data_table):
		bar_set = ax.bar(ind + (iteration * width),data_table[key],width,align="center")
		#print(data_table[key])
		iteration += 1
		# only use the sample name
		head, tail = os.path.split(key)
		# remove the mpileup_tab
		tail = tail.replace('_mpileup.tab','')
		legend_list.append(tail)

	plt.xlabel('Position in spacer')
	if modality=='BE4':
		plt.ylabel('C>T conversion frequency')
	else:
		plt.ylabel('A>G conversion frequency')

	# set the axis tick labels
	ax.legend(legend_list)

	# save to file
	plt.savefig(output_file)

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--mutation_summary_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-d','--diversity_summary_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--mutation_graph_output',required=True)
	parser.add_argument('-g','--diversity_graph_output',required=True)	
	opts = parser.parse_args(argv)
	if opts.modality=='ABE7.10' or opts.modality=='BE4':
		make_base_editing_plot(opts.input_file,opts.modality,opts.output_file)
	else:
		make_indel_rate_plot(opts.mutation_summary_file,opts.mutation_graph_output)
		make_diversity_plot(opts.diversity_summary_file,opts.diversity_graph_output)
 
if __name__ == '__main__':
	main(sys.argv[1:])
