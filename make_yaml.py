# script to make a yaml file from a tab-delimited sample file
# in the following order
# Sample name: in the format of 'Plate-Well-Gene'
# Reference file: reference sequence of the amplicon being sequenced
# sgRNA target site: target site for the CRISPR experiment
# Experiment type: HDR/KO
# if HDR: variants (semicolon delimited, no spaces), else NA
# Control sample: PCR product from an untreated sample (or appropriate control)




# output file: yaml file which can run in the pipeline

# input file for genotyping: sample name {tab} reference {SNP diffs ; delimited A1008T}

# import appropriate libraries
import yaml
import sys
import argparse
import re
import glob
import subprocess
from collections import defaultdict

def makeYAML (input_file,project_name,ngs_run,cell_type,output_file,directory_name):
	# dictionary to hold yaml
	yaml_dict = {}
	yaml_dict['samples'] = []
	yaml_dict['project'] = project_name
	yaml_dict['ngs_run'] = ngs_run
	yaml_dict['cell_type'] = cell_type
	ctrls_done = []

	for line in input_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		#print(line)

		# variables to store fields
		sample_name = parts[0]
		reference_file = parts[1]
		coordinates = parts[2]
		target_site = parts[3]
		expt_type = parts[4]
		variant_list = parts[5]
		control_sample = parts[6]

		# save to the yaml
		sample_dict = {}
		sample_dict[sample_name] = {}
		sample_dict[sample_name]['name'] = sample_name
		sample_dict[sample_name]['reference'] = reference_file
		sample_dict[sample_name]['coordinates'] = coordinates
		sample_dict[sample_name]['target_site'] = target_site
		sample_dict[sample_name]['expt_type'] = expt_type
		sample_dict[sample_name]['variants'] = variant_list
		sample_dict[sample_name]['control'] = control_sample

		# add to yaml dict
		yaml_dict['samples'].append(sample_dict)
		
		# go through and change the filenames
		path_to_check = directory_name + '/' + sample_name + '_*.fastq.gz'
		file_list = sorted(glob.glob(path_to_check))

		# and change each file
		mvCommand_R1 = 'mv ' + file_list[0] + ' ' + directory_name + '/' + sample_name + '_R1.fastq.gz'
		p = subprocess.Popen(mvCommand_R1,shell=True)
		p.communicate()
		mvCommand_R2 = 'mv ' + file_list[1] + ' ' + directory_name + '/' + sample_name + '_R2.fastq.gz'
		p = subprocess.Popen(mvCommand_R2,shell=True)
		p.communicate()

		# do this for controls
		if control_sample not in ctrls_done:
			path_to_check = directory_name + '/' + control_sample + '_*.fastq.gz'
			file_list = sorted(glob.glob(path_to_check))
			mvCommand_R1 = 'mv ' + file_list[0] + ' ' + directory_name + '/' + control_sample + '_R1.fastq.gz'
			p = subprocess.Popen(mvCommand_R1,shell=True)
			p.communicate()
			mvCommand_R2 = 'mv ' + file_list[1] + ' ' + directory_name + '/' + control_sample + '_R2.fastq.gz'
			p = subprocess.Popen(mvCommand_R2,shell=True)
			p.communicate()
			ctrls_done.append(control_sample)

	# write final yaml file
	with open(output_file, 'w') as yaml_file:
		yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

    # close everything
	input_file.close()
	yaml_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-p','--project_name',required=True)
	parser.add_argument('-n','--ngs_run',required=True)
	parser.add_argument('-c','--cell_type',required=True)
	parser.add_argument('-d','--directory_name',required=True)
	parser.add_argument('-o','--output_file',required=True)
	opts = parser.parse_args(argv)
	makeYAML(opts.input_file,opts.project_name,opts.ngs_run,opts.cell_type,opts.output_file,opts.directory_name)
 
if __name__ == '__main__':
	main(sys.argv[1:])