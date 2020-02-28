import os

configfile: "pan_config.json"

fastq_directory = "/mnt/storage/COMBINEDLABPUB/greatapes/fastq"

temp_directory = "temp/"

gatk_path = "/home/thwebste/Tools/GenomeAnalysisTK_37.jar"

bbduksh_path = "bbduk.sh"
bbmerge_sh_path = "bbmerge.sh"
bcftools_path = "bcftools"
bgzip_path = "bgzip"
bwa_path = "bwa"
fastq_dump_path = "fastq-dump"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
picard_path = "picard"
prefetch_path = "prefetch"
sambamba_path = "sambamba"
samtools_path = "samtools"
tabix_path = "tabix"

paired = [x for x in config["sras"] if x not in config["single_end"]]

# exclude_list = ["SRR740818", "SRR740831"]
exclude_list = []

fastq_prefixes_paired = [
	config[x]["fq1"][:-9] for x in paired] + [
		config[x]["fq2"][:-9] for x in paired]

fastq_prefixes_single = [
	config[x]["fq1"][:-9] for x in config["single_end"]]

trimmed_fastq_prefixes = [
	"{}_trimmed_read1".format(x) for x in paired] + [
	"{}_trimmed_read2".format(x) for x in paired] + [
	"{}_trimmed_single".format(x) for x in config["single_end"]]

sample_dictionary = {}

for k in config["samples"]:
	sample_dictionary[k] = [
		x for x in config["samples"][k] if x not in exclude_list]

paired_samples = {}
for k in sample_dictionary:
	paired_samples[k] = [
		x for x in sample_dictionary[k] if x not in config["single_end"]]

single_samples = {}
for k in sample_dictionary:
	single_samples[k] = [
		x for x in sample_dictionary[k] if x in config["single_end"]]


def gather_bam_input_for_merging(wildcards):
	if len(single_samples[wildcards.sample]) > 0:
		singles_list = expand(
			"processed_bams/{sample}.{genome}.sorted.single.bam",
			sample=single_samples[wildcards.sample], genome=wildcards.genome)
	else:
		singles_list = []
	paired_list = expand(
		"processed_bams/{sample}.{genome}.sorted.paired.bam",
		sample=paired_samples[wildcards.sample], genome=wildcards.genome)
	return singles_list + paired_list


def gather_bai_input_for_merging(wildcards):
	if len(single_samples[wildcards.sample]) > 0:
		singles_list = expand(
			"processed_bams/{sample}.{genome}.sorted.single.bam.bai",
			sample=single_samples[wildcards.sample], genome=wildcards.genome)
	else:
		singles_list = []
	paired_list = expand(
		"processed_bams/{sample}.{genome}.sorted.paired.bam.bai",
		sample=paired_samples[wildcards.sample], genome=wildcards.genome)
	return singles_list + paired_list


assemblies = ["pantro6"]

rule all:
	input:
		expand(
			"xyalign/reference/{assembly}.{ver}.fasta.fai",
			assembly=assemblies, ver=["XY", "XXonly"]),
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html",
		# expand(
		# 	"stats/{sample}.{genome}.sorted.mkdup.bam.stats",
		# 	sample=config["sample_names"], genome=assemblies),
		# expand(
		# 	"callable_sites/{sample}.{genome}.ONLYcallablesites.bed",
		# 	sample=config["sample_names"], genome=assemblies),
		# expand(
		# 	"vcf_combined/{population}.{genome}.combined.filtered_{type}.vcf.gz.tbi",
		# 	population=["allpan"],
		# 	genome=assemblies,
		# 	type=["allvariant"])

# rule prepare_reference:
# 	input:
# 		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
# 	output:
# 		new = "reference/{assembly}.fasta",
# 		fai = "reference/{assembly}.fasta.fai",
# 		amb = "reference/{assembly}.fasta.amb",
# 		dict = "reference/{assembly}.dict"
# 	params:
# 		samtools = samtools_path,
# 		bwa = bwa_path
# 	run:
# 		shell(
# 			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
# 		# faidx
# 		shell(
# 			"{params.samtools} faidx {output.new}")
# 		# .dict
# 		shell(
# 			"{params.samtools} dict -o {output.dict} {output.new}")
# 		# bwa
# 		shell(
# 			"{params.bwa} index {output.new}")

rule prefetch_sra:
	output:
		os.path.join(temp_directory, "{id}/{id}.sra")
	params:
		tool = prefetch_path,
		tmp_dir = temp_directory,
		use_id = "{id}"
	shell:
		"{params.tool} {params.use_id} -O {params.tmp_dir}"

rule fastq_dump_paired:
	input:
		sra = os.path.join(temp_directory, "{sample}/{sample}.sra")
	output:
		fq1 = os.path.join("paired_fastqs", "{sample}_1.fastq.gz"),
		fq2 = os.path.join("paired_fastqs", "{sample}_2.fastq.gz")
	params:
		output_dir = "paired_fastqs",
		fastq_dump = fastq_dump_path
	shell:
		"{params.fastq_dump} --outdir {params.output_dir} --gzip --readids --split-files {input.sra}"

rule fastq_dump_single:
	input:
		sra = os.path.join(temp_directory, "{sample}/{sample}.sra")
	output:
		fq1 = os.path.join("single_fastqs", "{sample}_1.fastq.gz")
	params:
		output_dir = "single_fastqs",
		fastq_dump = fastq_dump_path
	shell:
		"{params.fastq_dump} --outdir {params.output_dir} --gzip --readids --split-files {input.sra}"

rule fastqc_analysis_paired:
	input:
		os.path.join("paired_fastqs", "{fq_prefix}.fastq.gz")
	output:
		"fastqc_paired/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule fastqc_analysis_single:
	input:
		os.path.join("single_fastqs", "{fq_prefix}.fastq.gz")
	output:
		"fastqc_single/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		paired = expand("fastqc_paired/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes_paired),
		single = expand("fastqc_single/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes_single)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc fastqc_paired fastqc_single"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			"paired_fastqs", config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			"paired_fastqs", config[wildcards.sample]["fq2"])
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

rule trim_adapters_single_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			"single_fastqs", config[wildcards.sample]["fq1"])
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_single.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in={input.fq1} "
		"out={output.out_fq1} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		"trimmed_fastqs/{fq_prefix}.fastq.gz"
	output:
		"trimmed_fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o trimmed_fastqc {input}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"trimmed_fastqc/{fq_prefix}_fastqc.html",
			fq_prefix=trimmed_fastq_prefixes)
	output:
		"multiqc_trimmed/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed trimmed_fastqc"

rule get_annotation:
	output:
		"reference/{genome}.fa"
	params:
		web_address = lambda wildcards: config["genome_paths"][wildcards.genome],
		initial_output = "reference/{genome}.fa.gz"
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("gunzip {params.initial_output}")

rule xyalign_create_references:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		xx = "xyalign/reference/{assembly}.XXonly.fa",
		xy = "xyalign/reference/{assembly}.XY.fasta"
	conda:
		"envs/xyalign.yml"
	params:
		gen_assembly = "{assembly}",
		output_dir = "xyalign",
		x = lambda wildcards: config["chrx"][wildcards.assembly],
		y = lambda wildcards: config["chry"][wildcards.assembly]
	shell:
		"xyalign --PREPARE_REFERENCE --ref {input} --bam null "
		"--xx_ref_out {output.xx} --xy_ref_out {output.xy} "
		"--output_dir {params.output_dir} --x_chromosome {params.x} "
		"--y_chromosome {params.y}"

rule prepare_reference_males:
	input:
		"xyalign/reference/{assembly}.XY.fasta"
	output:
		fai = "xyalign/reference/{assembly}.XY.fasta.fai",
		amb = "xyalign/reference/{assembly}.XY.fasta.amb",
		dict = "xyalign/reference/{assembly}.XY.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

rule prepare_reference_females:
	input:
		"xyalign/reference/{assembly}.XXonly.fa"
	output:
		fai = "xyalign/reference/{assembly}.XXonly.fa.fai",
		amb = "xyalign/reference/{assembly}.XXonly.fa.amb",
		dict = "xyalign/reference/{assembly}.XXonly.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

rule map_and_process_trimmed_paired_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		fai_xx = "xyalign/reference/{assembly}.XXonly.fa.fai",
		ref_xx = "xyalign/reference/{assembly}.XXonly.fa",
		fai_xy = "xyalign/reference/{assembly}.XY.fasta.fai",
		ref_xy = "xyalign/reference/{assembly}.XY.fasta"
	output:
		"processed_bams/{sample}.{assembly}.sorted.paired.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path
	threads: 4
	run:
		if wildcards.sample in males:
			shell(
				" {params.bwa} mem -t {threads} -R "
			 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
				"{input.ref_xy} {input.fq1} {input.fq2}"
				"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
				"-O bam -o {output}")
		else:
			shell(
				" {params.bwa} mem -t {threads} -R "
			 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
				"{input.ref_xx} {input.fq1} {input.fq2}"
				"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
				"-O bam -o {output}")

rule map_and_process_trimmed_single_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_single.fastq.gz",
		fai_xx = "xyalign/reference/{assembly}.XXonly.fa.fai",
		ref_xx = "xyalign/reference/{assembly}.XXonly.fa",
		fai_xy = "xyalign/reference/{assembly}.XY.fasta.fai",
		ref_xy = "xyalign/reference/{assembly}.XY.fasta"
	output:
		"processed_bams/{sample}.{assembly}.sorted.single.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path
	threads: 4
	run:
		if wildcards.sample in males:
			shell(
				" {params.bwa} mem -t {threads} -R "
			 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
				"{input.ref_xy} {input.fq1} "
				"| {params.samtools} view -b - | {params.samtools} sort "
				"-O bam -o {output}")
		else:
			shell(
				" {params.bwa} mem -t {threads} -R "
			 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
				"{input.ref_xx} {input.fq1} "
				"| {params.samtools} view -b - | {params.samtools} sort "
				"-O bam -o {output}")

rule index_paired_bam:
	input:
		"processed_bams/{sample}.{assembly}.sorted.paired.bam"
	output:
		"processed_bams/{sample}.{assembly}.sorted.paired.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule index_single_bam:
	input:
		"processed_bams/{sample}.{assembly}.sorted.single.bam"
	output:
		"processed_bams/{sample}.{assembly}.sorted.single.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule merge_bams:
	input:
		bams = gather_bam_input_for_merging,
		bais = gather_bai_input_for_merging
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

rule generate_callable_sites:
	input:
		ref = "reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"callable_sites/{sample}.{genome}.callablesites"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		summary = "stats/{sample}.{genome}.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} "
		"-jar {params.gatk} -T CallableLoci -R {input.ref} "
		"-I {input.bam} --minDepth 10 --minMappingQuality 30 "
		"--summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.{genome}.callablesites"
	output:
		"callable_sites/{sample}.{genome}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"

# rule gatk_gvcf:
# 	input:
# 		ref = "reference/{genome}.fasta",
# 		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
# 		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
# 	output:
# 		"vcf/{sample}.{genome}.g.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path
# 	threads:
# 		4
# 	shell:
# 		"java -Xmx15g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} "
# 		"-T HaplotypeCaller -R {input.ref} -I {input.bam} "
# 		"-contamination 0.05 --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_per_chrom:
	input:
		ref = "reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"vcf/{sample}.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		4
	shell:
		"java -Xmx15g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} "
		"-T HaplotypeCaller -R {input.ref} -I {input.bam} -L {wildcards.chrom} "
		"-contamination 0.05 --emitRefConfidence GVCF -o {output}"

rule genotype_gvcfs_per_chrom:
	input:
		ref = "reference/{genome}.fasta",
		gvcfs = lambda wildcards: expand(
			"vcf/{sample}.{genome}.{chrom}.g.vcf.gz",
			sample=config["subspecies"][wildcards.population],
			genome=[wildcards.genome], chrom=[wildcards.chrom])
	output:
		v = "vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"java -Xmx15g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} "
			"-T GenotypeGVCFs -R {input.ref} {variant_files} "
			"-o {output.v}")

rule filter_vcfs_polymorphic_only:
	input:
		"vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
	output:
		"vcf_joint/{population}.{genome}.{chrom}.filtered_polymorphic.vcf.gz"
	params:
		bgzip = bgzip_path,
		bcftools = bcftools_path
	shell:
		"{params.bcftools} filter -i "
		"'QUAL >= 30 && MQ >= 30 && AF > 0 && AF < 1.0 && QD > 2' {input} | "
		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
		"{params.bgzip} > {output}"

rule index_filtered_vcf_polymorphic:
	input:
		"vcf_joint/{population}.{genome}.{chrom}.filtered_polymorphic.vcf.gz"
	output:
		"vcf_joint/{population}.{genome}.{chrom}.filtered_polymorphic.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule filter_vcfs_allvariant:
	input:
		"vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
	output:
		"vcf_joint/{population}.{genome}.{chrom}.filtered_allvariant.vcf.gz"
	params:
		bgzip = bgzip_path,
		bcftools = bcftools_path
	shell:
		"{params.bcftools} filter -i "
		"'QUAL >= 30 && MQ >= 30 && QD > 2' {input} | "
		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
		"{params.bgzip} > {output}"

rule index_filtered_vcf_allvariant:
	input:
		"vcf_joint/{population}.{genome}.{chrom}.filtered_allvariant.vcf.gz"
	output:
		"vcf_joint/{population}.{genome}.{chrom}.filtered_allvariant.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule concatenate_split_vcfs:
	input:
		vcf = lambda wildcards: expand(
			"vcf_joint/{pop}.{gen}.{chrom}.filtered_{ty}.vcf.gz",
			pop=wildcards.population,
			gen=wildcards.genome,
			ty=wildcards.type,
			chrom=config["chromosomes_to_analyze"][wildcards.genome]),
		idx = lambda wildcards: expand(
			"vcf_joint/{pop}.{gen}.{chrom}.filtered_{ty}.vcf.gz.tbi",
			pop=wildcards.population,
			gen=wildcards.genome,
			ty=wildcards.type,
			chrom=config["chromosomes_to_analyze"][wildcards.genome])
	output:
		"vcf_combined/{population}.{genome}.combined.filtered_{type}.vcf.gz"
	params:
		bcftools = bcftools_path
	shell:
		"{params.bcftools} concat -O z -o {output} {input.vcf}"

rule index_concatenated_vcf:
	input:
		"vcf_combined/{population}.{gen}.combined.filtered_{type}.vcf.gz"
	output:
		"vcf_combined/{population}.{gen}.combined.filtered_{type}.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

# rule filter_vcfs:
# 	input:
# 		"vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
# 	output:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered.vcf.gz"
# 	params:
# 		bgzip = bgzip_path,
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} filter -i "
# 		"'QUAL >= 30 && MQ >= 30 && AF > 0 && AF < 1.0 && QD > 2 & "
# 		"FMT/DP >= 10 & FMT/GQ >= 30' {input} | "
# 		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
# 		"{params.bgzip} > {output}"
#
# rule index_filtered_vcf:
# 	input:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered.vcf.gz"
# 	output:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"
