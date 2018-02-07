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
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html"

rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
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

rule trim_adapters_paired_bbduk:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}_read1.fastq.gz"),
		fq2 = os.path.join(fastq_directory, "{sample}_read2.fastq.gz")
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "trimmed_fastqc/{sample}_trimmed_read1_fastqc.html",
		html2 = "trimmed_fastqc/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o trimmed_fastqc {input.fq1} {input.fq2}"

rule multiqc_analysis:
	input:
		expand(
			"trimmed_fastqc/{sample}_trimmed_{read}_fastqc.html",
			sample=config["sras"], read=["read1", "read2"])
	output:
		"multiqc_trimmed/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc_trimmed trimmed_fastqc"
