# script to make a yaml file from a tab-delimited sample file
# updated 09/18/18 to compute control subtracted rates

# input file: [tab-delimited - sample name {tab} reference file {tab} target site {tab} control_sample
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

def makeYAML (input_file,project_name,ngs_run,cell_type,output_file,directory_name,pipeline_type):
	# dictionary to hold yaml
	yaml_dict = {}
	yaml_dict['samples'] = []
	yaml_dict['project'] = project_name
	yaml_dict['ngs_run'] = ngs_run
	yaml_dict['cell_type'] = cell_type

	for line in input_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')

		# variables to store fields
		sample_name = parts[0]
		reference_file = parts[1]
		variant_list = parts[2]

		# save to the yaml
		sample_dict = {}
		sample_dict[sample_name] = {}
		sample_dict[sample_name]['name'] = sample_name
		sample_dict[sample_name]['reference'] = reference_file
		sample_dict[sample_name]['variants'] = variant_list

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