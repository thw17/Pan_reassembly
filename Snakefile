import os

configfile: "pan_config.json"

fastq_directory = "/mnt/storage/COMBINEDLABPUB/greatapes/fastq"

temp_directory = "temp/"

bbduksh_path = "bbduk.sh"
bbmerge_sh_path = "bbmerge.sh"
bwa_path = "bwa"
fastqc_path = "fastqc"
gatk_path = "gatk-launch"
multiqc_path = "multiqc"
picard_path = "picard"
sambamba_path = "sambamba"
samtools_path = "samtools"

fastq_prefixes = [
	config[x]["fq1"][:-9] for x in config["sras"]] + [
		config[x]["fq2"][:-9] for x in config["sras"]]

rule all:
	input:
		expand("reference/{assembly}.fasta.fai", assembly=["pantro4"]),
		expand("adapters/{sample}.adapters.fa", sample=config["sras"]),
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html",
		expand(
			"processed_bams/{sample}.{genome}.sorted.merged.bam.bai",
			sample=config["sample_names"], genome=["pantro4"]),
		expand(
			"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai",
			sample=config["sample_names"], genome=["pantro4"]),
		expand(
			"stats/{sample}.{genome}.sorted.mkdup.bam.stats",
			sample=config["sample_names"], genome=["pantro4"]),
		expand(
			"vcf/{sample}.{genome}.g.vcf.gz",
			sample=config["sample_names"], genome=["pantro4"]),

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
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
	output:
		"adapters/{sample}.adapters.fa"
	params:
		bbmerge_sh = bbmerge_sh_path
	shell:
		"{params.bbmerge_sh} in1={input.fq1} in2={input.fq2} outa={output} reads=1m"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
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

rule multiqc_analysis_trimmed:
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

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		new = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"processed_bams/{sample}.{assembly}.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	shell:
		" {params.bwa} mem -t {params.threads} -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.new} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule merge_bams:
	input:
		bams = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.sorted.bam",
			sample=config["samples"][wildcards.sample], genome=wildcards.genome),
		bais = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.sorted.bam.bai",
			sample=config["samples"][wildcards.sample], genome=wildcards.genome)
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	threads: 4
	params:
		sambamba = sambamba_path,
		threads = 4
	shell:
		"{params.sambamba} merge -t {params.threads} {output} {input.bams}"

rule index_merged_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

# rule sambamba_mark_dups_only_analysis_chroms:
# 	input:
# 		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
# 		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
# 	output:
# 		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
# 	threads: 4
# 	params:
# 		sambamba = sambamba_path,
# 		samtools = samtools_path,
# 		regions = lambda wildcards: config[
# 			"chromosomes_to_analyze"][wildcards.genome],
# 		threads = 4
# 	shell:
# 		"{params.samtools} view -b {input.bam} {params.regions} | "
# 		"{params.sambamba} markdup -t {params.threads} /dev/stdin {output}"

rule picard_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		metrics = "stats/{sample}.{genome}.picard_mkdup_metrics.txt"
	threads: 4
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule gatk_gvcf:
	input:
		ref = "reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"vcf/{sample}.{genome}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		4
	shell:
		"{params.gatk} --java-options '-Xmx15g -Djava.io.tmpdir={params.temp_dir}' HaplotypeCaller -R {input.ref} -I {input.bam} -contamination 0.05 --emit-ref-confidence GVCF -o {output}"
