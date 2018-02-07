import os

configfile: "pan_config.json"

fastq_directory = "/mnt/storage/COMBINEDLABPUB/greatapes/fastq"

temp_directory = "temp/"

bbmerge_sh_path = "bbmerge.sh"
bwa_path = "bwa"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
samtools_path = "samtools"

fastq_prefixes = [
	config[x]["fq1"][:-9] for x in config["sras"]] + [
		config[x]["fq2"][:-9] for x in config["sras"]]

rule all:
	input:
		expand("reference/{assembly}.fasta.fai", assembly=["pantro4"]),
		expand("adapters/{sample}.adapters.fa", sample=config["sras"]),
		"multiqc/multiqc_report"

rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "new_reference/{assembly}.fasta",
		fai = "new_reference/{assembly}.fasta.fai",
		amb = "new_reference/{assembly}.fasta.amb",
		dict = "new_reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		shell(
			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell(
			"{params.samtools} faidx {output.new}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {output.new}")
		# bwa
		shell(
			"{params.bwa} index {output.new}")

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc fastqc"

rule adapter_discovery:
	input:
		fq1 = lambda wildcards: os.path.join(fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(fastq_directory, config[wildcards.sample]["fq2"])
	output:
		"adapters/{sample}.adapters.fa"
	params:
		bbmerge_sh = bbmerge_sh_path
	shell:
		"{params.bbmerge_sh} in1={input.fq1} in2={input.fq2} outa={output} reads=1m"
