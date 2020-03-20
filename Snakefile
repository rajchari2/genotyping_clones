# test

import snakemake
import os
from collections import defaultdict

# get all of the samples
samples = config["samples"]
project_name = config["project"]
cell_type = config["cell_type"]
ngs_run = config["ngs_run"]
sample_list = []
reference_list = []
control_list = []
sample_to_reference = defaultdict(str)
sample_to_vars = defaultdict(str)
sample_to_control = defaultdict(str)
sample_to_target_site = defaultdict(str)
sample_to_expt_type = defaultdict(str)

for sample in samples:
	for sample_name in sample:
		sample_list.append(sample_name)
		ref_file = sample[sample_name]['reference']
		if ref_file not in reference_list:
			reference_list.append(ref_file)
		sample_to_reference[sample_name] = 'db/' + ref_file
		sample_to_vars[sample_name] = sample[sample_name]['variants']
		sample_to_expt_type[sample_name] = sample[sample_name]['expt_type']
		sample_to_control[sample_name] = sample[sample_name]['control']
		sample_to_target_site[sample_name] = sample[sample_name]['target_site']
		control_sample = sample[sample_name]['control']
		if control_sample not in control_list:
			control_list.append(control_sample)
			sample_to_reference[control_sample] = 'db/' + ref_file

# functions to get variables
def get_reference(wildcards):
	return sample_to_reference[wildcards.sample]

def get_target_site(wildcards):
	return sample_to_target_site[wildcards.sample]

def get_control_sample(wildcards):
	return sample_to_control[wildcards.sample]

def get_variants(wildcards):
	return sample_to_vars[wildcards.sample]

def get_expt_type(wildcards):
	return sample_to_expt_type[wildcards.sample]

# rules for genotyping
rule bwa_mem:
	input:
		reference_file = get_reference,
		r1 = 'data_files/{sample}_R1.fastq.gz',
		r2 = 'data_files/{sample}_R2.fastq.gz',
	output:
		sam_file = 'processed/{sample}_bwamem.sam'
	shell:
		'bwa index {input.reference_file} && bwa mem {input.reference_file} {input.r1} {input.r2} > {output.sam_file}'

rule samtools_view:
	input:
		sam = rules.bwa_mem.output.sam_file,
		reference_file = get_reference
	output:
		bam = 'processed/{sample}_bwamem.bam'
	shell:
		'samtools view -bt {input.reference_file} -o {output.bam} {input.sam}'

rule samtools_sort:
	input:
		bam = rules.samtools_view.output.bam
	output:
		sorted_bam = 'processed/{sample}_bwamem_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.bam}'

rule bam_index:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam
	output:
		bai = 'processed/{sample}_bwamem_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

rule count_variants:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam,
		bam_index = rules.bam_index.output.bai
	params:
		variant_list = get_variants,
		expt_type = get_expt_type,
		control = get_control_sample,
		reference_file = get_reference,
		target_site = get_target_site
	output:
		var_file = 'output/{sample}_variant_summary.tab'
	shell:
		'python resources/genotype_sample.py -i {input.sorted_bam} -c {params.control} -e {params.expt_type} -v {params.variant_list} -o {output.var_file} -r {params.reference_file} -t {params.target_site}'

rule generate_toml_file:
	input:
		variant_file_list = expand(rules.count_variants.output.var_file,sample=sample_list)
	output:
		toml_file = 'final_output/' + project_name + '_' + ngs_run + '_genotyping.toml'
	shell:
		'python resources/aggregate_genotypes.py -i {input.variant_file_list} -o {output.toml_file}'

rule visualize_plate:
	input:
		toml_file = rules.generate_toml_file.output.toml_file
	params:
		color = 'Greys'
	output:
		plate_map = 'final_output/' + project_name + '_' + ngs_run + '_clonotyping_plot.png'
	shell:
		'bio96 -o {output.plate_map} {input.toml_file} -c {params.color}'

rule build_controls:
	input:
		pileup_list = expand(rules.samtools_sort.output.sorted_bam,sample=control_list),
		bais = expand(rules.bam_index.output.bai,sample=control_list)
	output:
		controls_created = project_name + '_' + ngs_run + '_Controls.tab'
	shell:
		'touch {output.controls_created}'

rule genotype_all:
	input:
		plate_map = rules.visualize_plate.output.plate_map


