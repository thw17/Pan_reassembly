import os

configfile: "pan_config.json"

temp_directory = "temp/"

# Runtime values
very_short = "6:00:00"
medium = "12:00:00"
day = "24:00:00"
long = "48:00:00"

# Tool paths
bbduksh_path = "bbduk.sh"
bbmerge_sh_path = "bbmerge.sh"
bcftools_path = "bcftools"
bgzip_path = "bgzip"
bwa_path = "bwa"
fastq_dump_path = "fastq-dump"
fastqc_path = "fastqc"
gatk_path = "gatk"
mosdepth_path = "mosdepth"
multiqc_path = "multiqc"
picard_path = "picard"
prefetch_path = "prefetch"
rename_sh_path = "rename.sh"
sambamba_path = "/uufs/chpc.utah.edu/common/home/u6023206/Programs/sambamba-0.7.1"
samtools_path = "samtools"
tabix_path = "tabix"

paired = [x for x in config["sras"] if x not in config["single_end"]]

males = [x for x in config["sample_names"] if x not in config["females"]]
females = config["females"]

exclude_list = ["SRR740818", "SRR740831"]

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
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html",
		# expand(
		# 	"xyalign/reference/{assembly}.{ver}.fa.fai",
		# 	assembly=assemblies, ver=["XY", "XXonly"]),
		expand(
			"stats/{sample}.{genome}.sorted.mkdup.bam.{tool}.stats",
			sample=config["sample_names"], genome=assemblies, tool=["picard", "sambamba"]),
		# expand(
		# 	"vcf_genotyped/pantro6.{chrom}.gatk.called.raw.vcf.gz",
		# 	chrom=config["chromosomes_to_analyze"]["pantro6"]),
		expand(
			"contamination_filter_vcf_genotyped/pantro6.{chrom}.gatk.called.raw.vcf.gz",
			chrom=config["chromosomes_to_analyze"]["pantro6"]),
		expand(
			"mosdepth_results/{sample}.{genome}.total.mosdepth.summary.txt",
			sample=config["sample_names"], genome=assemblies)
		# expand(
		# 	"callable_sites/{sample}.{genome}.ONLYcallablesites.bed",
		# 	sample=config["sample_names"], genome=assemblies),
		# expand(
		# 	"vcf_combined/{population}.{genome}.combined.filtered_{type}.vcf.gz.tbi",
		# 	population=["allpan"],
		# 	genome=assemblies,
		# 	type=["allvariant"])

rule prefetch_sra:
	output:
		os.path.join(temp_directory, "{id}/{id}.sra")
	params:
		tool = prefetch_path,
		tmp_dir = temp_directory,
		use_id = "{id}",
		threads = 1,
		mem = 4,
		t = day
	shell:
		"{params.tool} {params.use_id} -O {params.tmp_dir} --max-size 100GB"

rule fastq_dump_paired:
	input:
		sra = os.path.join(temp_directory, "{sample}/{sample}.sra")
	output:
		fq1 = os.path.join("paired_fastqs", "{sample}_1.fastq.gz"),
		fq2 = os.path.join("paired_fastqs", "{sample}_2.fastq.gz")
	params:
		output_dir = "paired_fastqs",
		fastq_dump = fastq_dump_path,
		threads = 1,
		mem = 4,
		t = day
	shell:
		"{params.fastq_dump} --outdir {params.output_dir} --gzip --readids --split-files {input.sra}"

rule fastq_dump_single:
	input:
		sra = os.path.join(temp_directory, "{sample}/{sample}.sra")
	output:
		fq1 = os.path.join("single_fastqs", "{sample}_1.fastq.gz")
	params:
		output_dir = "single_fastqs",
		fastq_dump = fastq_dump_path,
		threads = 1,
		mem = 4,
		t = day
	shell:
		"{params.fastq_dump} --outdir {params.output_dir} --gzip --readids --split-files {input.sra}"

rule fastqc_analysis_paired:
	input:
		os.path.join("paired_fastqs", "{fq_prefix}.fastq.gz")
	output:
		"fastqc_paired/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc_paired {input}"

rule fastqc_analysis_single:
	input:
		os.path.join("single_fastqs", "{fq_prefix}.fastq.gz")
	output:
		"fastqc_single/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc_single {input}"

rule multiqc_analysis:
	input:
		paired = expand("fastqc_paired/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes_paired),
		single = expand("fastqc_single/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes_single)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		bbduksh = bbduksh_path,
		threads = 2,
		mem = 8,
		t = very_short
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
		bbduksh = bbduksh_path,
		threads = 2,
		mem = 8,
		t = very_short
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
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o trimmed_fastqc {input}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"trimmed_fastqc/{fq_prefix}_fastqc.html",
			fq_prefix=trimmed_fastq_prefixes)
	output:
		"multiqc_trimmed/multiqc_report.html",
	params:
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed trimmed_fastqc"

rule fix_read_IDs_for_paired_fastqs_from_SRA_paired:
	# The fastq-dump created issues with read ID names for paired files so that
	# they give bwa issues.  This rule will go through and rename them so that
	# they're compatible with bwa
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		out1 = "renamed_fastqs/{sample}_trimmed_fixed_1.fastq.gz",
		out2 = "renamed_fastqs/{sample}_trimmed_fixed_2.fastq.gz"
	params:
		rename_sh = rename_sh_path,
		read_name = "{sample}",
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.rename_sh} in={input.fq1} in2={input.fq2} out={output.out1} out2={output.out2} prefix={params.read_name}"

rule fix_read_IDs_for_paired_fastqs_from_SRA_single:
	# The fastq-dump created issues with read ID names for paired files so that
	# they give bwa issues.  This rule will go through and rename them so that
	# they're compatible with bwa
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_single.fastq.gz"
	output:
		out1 = "renamed_fastqs/{sample}_trimmed_fixed_single.fastq.gz"
	params:
		rename_sh = rename_sh_path,
		read_name = "{sample}",
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.rename_sh} in={input.fq1} out={output.out1} prefix={params.read_name}"

rule get_reference:
	output:
		"reference/{genome}.fa"
	params:
		web_address = lambda wildcards: config["genome_paths"][wildcards.genome],
		initial_output = "reference/{genome}.fa.gz",
		threads = 1,
		mem = 4,
		t = very_short
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("gunzip {params.initial_output}")

rule xyalign_create_references:
	input:
		ref = "reference/{assembly}.fa"
	output:
		xx = "xyalign/reference/{assembly}.XXonly.fa",
		xy = "xyalign/reference/{assembly}.XY.fa"
	conda:
		"envs/xyalign.yml"
	params:
		gen_assembly = "{assembly}",
		output_dir = "xyalign",
		x = lambda wildcards: config["chrx"][wildcards.assembly],
		y = lambda wildcards: config["chry"][wildcards.assembly],
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"xyalign --PREPARE_REFERENCE --ref {input} --bam null "
		"--xx_ref_out {output.xx} --xy_ref_out {output.xy} "
		"--output_dir {params.output_dir} --x_chromosome {params.x} "
		"--y_chromosome {params.y}"

rule prepare_reference_males:
	input:
		"xyalign/reference/{assembly}.XY.fa"
	output:
		fai = "xyalign/reference/{assembly}.XY.fa.fai",
		amb = "xyalign/reference/{assembly}.XY.fa.amb",
		dict = "xyalign/reference/{assembly}.XY.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path,
		threads = 4,
		mem = 16,
		t = medium
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
		bwa = bwa_path,
		threads = 4,
		mem = 16,
		t = medium
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

rule map_and_process_trimmed_paired_reads:
	input:
		fq1 = "renamed_fastqs/{sample}_trimmed_fixed_1.fastq.gz",
		fq2 = "renamed_fastqs/{sample}_trimmed_fixed_2.fastq.gz",
		fai_xx = "xyalign/reference/{assembly}.XXonly.fa.fai",
		ref_xx = "xyalign/reference/{assembly}.XXonly.fa",
		fai_xy = "xyalign/reference/{assembly}.XY.fa.fai",
		ref_xy = "xyalign/reference/{assembly}.XY.fa"
	output:
		"processed_bams/{sample}.{assembly}.sorted.paired.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
	threads: 4
	run:
		if wildcards.sample in males:
			shell(
				" {params.bwa} mem -t {threads} -R "
			 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
				"{input.ref_xy} {input.fq1} {input.fq2}"
				"| {params.samtools} fixmate -m -O bam - - | {params.samtools} sort "
				"-O bam -o {output}")
		else:
			shell(
				" {params.bwa} mem -t {threads} -R "
			 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
				"{input.ref_xx} {input.fq1} {input.fq2}"
				"| {params.samtools} fixmate -m -O bam - - | {params.samtools} sort "
				"-O bam -o {output}")

rule map_and_process_trimmed_single_reads:
	input:
		fq1 = "renamed_fastqs/{sample}_trimmed_fixed_single.fastq.gz",
		fai_xx = "xyalign/reference/{assembly}.XXonly.fa.fai",
		ref_xx = "xyalign/reference/{assembly}.XXonly.fa",
		fai_xy = "xyalign/reference/{assembly}.XY.fa.fai",
		ref_xy = "xyalign/reference/{assembly}.XY.fa"
	output:
		"processed_bams/{sample}.{assembly}.sorted.single.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
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
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule index_single_bam:
	input:
		"processed_bams/{sample}.{assembly}.sorted.single.bam"
	output:
		"processed_bams/{sample}.{assembly}.sorted.single.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

# rule merge_bams:
# 	input:
# 		bams = gather_bam_input_for_merging,
# 		bais = gather_bai_input_for_merging
# 	output:
# 		"processed_bams/{sample}.{genome}.sorted.merged.bam"
# 	threads: 4
# 	params:
# 		sambamba = sambamba_path,
# 		threads = 4,
# 		mem = 16,
# 		t = long
# 	shell:
# 		"{params.sambamba} merge -t {params.threads} {output} {input.bams}"

rule merge_bams:
	input:
		bams = gather_bam_input_for_merging,
		bais = gather_bai_input_for_merging
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	threads: 4
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.samtools} merge {output} {input.bams}"

rule index_merged_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		metrics = "stats/{sample}.{genome}.picard_mkdup_metrics.txt"
	threads: 8
	params:
		picard = picard_path,
		threads = 8,
		mem = 24,
		t = long,
		tmp_dir = temp_directory
	shell:
		"{params.picard} -Xmx14g -Djava.io.tmpdir={params.tmp_dir} "
		"MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics} ASSUME_SORT_ORDER=coordinate"

rule index_picard_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats_picard:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.picard.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule mosdepth_total:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"mosdepth_results/{sample}.{genome}.total.mosdepth.summary.txt"
	params:
		mosdepth = mosdepth_path,
		prefix = "mosdepth_results/{sample}.{genome}.total",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.mosdepth} --fast-mode -F 1024 {params.prefix} {input.bam}"

rule samtools_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.samtools.mkdup.bam",
		metrics = "stats/{sample}.{genome}.samtools_mkdup_metrics.txt"
	threads: 8
	params:
		samtools = samtools_path,
		threads = 8,
		mem = 24,
		t = long,
		tmp_dir = temp_directory
	shell:
		"{params.samtools} markdup -s -f {output.metrics} --mode s - {output.bam}"

rule index_samtools_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.samtools.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.samtools.mkdup.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats_samtools:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.samtools.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.samtools.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.samtools.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule sambamba_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.sambamba.mkdup.bam"
	threads: 8
	params:
		sambamba = sambamba_path,
		threads = 8,
		mem = 32,
		t = long,
		tmp_dir = temp_directory
	shell:
		"{params.sambamba} markdup -t 6 --tmpdir={params.tmp_dir} {input.bam} {output.bam}"

rule index_sambamba_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.sambamba.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.sambamba.mkdup.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats_sambamba:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.sambamba.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.sambamba.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.sambamba.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule gatk_gvcf_per_chrom:
	input:
		ref = "xyalign/reference/{genome}.XY.fa",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"gvcf/{sample}.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		chr = "{chrom}",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chr} """
		"""-ERC GVCF --do-not-run-physical-phasing -O {output}"""

rule combine_gvcfs_per_chrom:
	input:
		ref = "xyalign/reference/{genome}.XY.fa",
		gvcfs = lambda wildcards: expand(
			"gvcf/{sample}.{genome}.{chrom}.g.vcf.gz",
			sample=config["sample_names"],
			genome=[wildcards.genome], chrom=[wildcards.chrom])
	output:
		v = "gvcf_combined/combined.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = long
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output.v}""")

rule gatk_genotypegvcf_per_chrom:
	input:
		ref = "xyalign/reference/{genome}.XY.fa",
		gvcf = "gvcf_combined/combined.{genome}.{chrom}.g.vcf.gz"
	output:
		"vcf_genotyped/{genome}.{chrom}.gatk.called.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""GenotypeGVCFs --include-non-variant-sites -R {input.ref} -V {input.gvcf} -O {output}"""


rule contamination_filter_gatk_gvcf_per_chrom:
	input:
		ref = "xyalign/reference/{genome}.XY.fa",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"contamination_filter_gvcf/{sample}.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		chr = "{chrom}",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chr} """
		"""-ERC GVCF --do-not-run-physical-phasing --contamination-fraction-to-filter 0.1 -O {output}"""

rule contamination_filter_combine_gvcfs_per_chrom:
	input:
		ref = "xyalign/reference/{genome}.XY.fa",
		gvcfs = lambda wildcards: expand(
			"contamination_filter_gvcf/{sample}.{genome}.{chrom}.g.vcf.gz",
			sample=config["sample_names"],
			genome=[wildcards.genome], chrom=[wildcards.chrom])
	output:
		v = "contamination_filter_gvcf_combined/combined.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = long
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output.v}""")

rule contamination_filter_gatk_genotypegvcf_per_chrom:
	input:
		ref = "xyalign/reference/{genome}.XY.fa",
		gvcf = "contamination_filter_gvcf_combined/combined.{genome}.{chrom}.g.vcf.gz"
	output:
		"contamination_filter_vcf_genotyped/{genome}.{chrom}.gatk.called.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""GenotypeGVCFs --include-non-variant-sites -R {input.ref} -V {input.gvcf} -O {output}"""

#
# rule generate_callable_sites:
# 	input:
# 		ref = "reference/{genome}.fasta",
# 		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
# 		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
# 	output:
# 		"callable_sites/{sample}.{genome}.callablesites"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path,
# 		summary = "stats/{sample}.{genome}.callable.summary"
# 	shell:
# 		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} "
# 		"-jar {params.gatk} -T CallableLoci -R {input.ref} "
# 		"-I {input.bam} --minDepth 10 --minMappingQuality 30 "
# 		"--summary {params.summary} -o {output}"
#
# rule extract_callable_sites:
# 	input:
# 		"callable_sites/{sample}.{genome}.callablesites"
# 	output:
# 		"callable_sites/{sample}.{genome}.ONLYcallablesites.bed"
# 	shell:
# 		"sed -e '/CALLABLE/!d' {input} > {output}"
#
# rule gatk_gvcf_per_chrom:
# 	input:
# 		ref = "reference/{genome}.fasta",
# 		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
# 		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
# 	output:
# 		"vcf/{sample}.{genome}.{chrom}.g.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path
# 	threads:
# 		4
# 	shell:
# 		"java -Xmx15g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} "
# 		"-T HaplotypeCaller -R {input.ref} -I {input.bam} -L {wildcards.chrom} "
# 		"-contamination 0.05 --emitRefConfidence GVCF -o {output}"
#
# rule genotype_gvcfs_per_chrom:
# 	input:
# 		ref = "reference/{genome}.fasta",
# 		gvcfs = lambda wildcards: expand(
# 			"vcf/{sample}.{genome}.{chrom}.g.vcf.gz",
# 			sample=config["subspecies"][wildcards.population],
# 			genome=[wildcards.genome], chrom=[wildcards.chrom])
# 	output:
# 		v = "vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path
# 	threads:
# 		4
# 	run:
# 		variant_files = []
# 		for i in input.gvcfs:
# 			variant_files.append("--variant " + i)
# 		variant_files = " ".join(variant_files)
# 		shell(
# 			"java -Xmx15g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} "
# 			"-T GenotypeGVCFs -R {input.ref} {variant_files} "
# 			"-o {output.v}")
#
# rule filter_vcfs_polymorphic_only:
# 	input:
# 		"vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
# 	output:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered_polymorphic.vcf.gz"
# 	params:
# 		bgzip = bgzip_path,
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} filter -i "
# 		"'QUAL >= 30 && MQ >= 30 && AF > 0 && AF < 1.0 && QD > 2' {input} | "
# 		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
# 		"{params.bgzip} > {output}"
#
# rule index_filtered_vcf_polymorphic:
# 	input:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered_polymorphic.vcf.gz"
# 	output:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered_polymorphic.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"
#
# rule filter_vcfs_allvariant:
# 	input:
# 		"vcf_joint/{population}.{genome}.{chrom}.vcf.gz"
# 	output:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered_allvariant.vcf.gz"
# 	params:
# 		bgzip = bgzip_path,
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} filter -i "
# 		"'QUAL >= 30 && MQ >= 30 && QD > 2' {input} | "
# 		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
# 		"{params.bgzip} > {output}"
#
# rule index_filtered_vcf_allvariant:
# 	input:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered_allvariant.vcf.gz"
# 	output:
# 		"vcf_joint/{population}.{genome}.{chrom}.filtered_allvariant.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"
#
# rule concatenate_split_vcfs:
# 	input:
# 		vcf = lambda wildcards: expand(
# 			"vcf_joint/{pop}.{gen}.{chrom}.filtered_{ty}.vcf.gz",
# 			pop=wildcards.population,
# 			gen=wildcards.genome,
# 			ty=wildcards.type,
# 			chrom=config["chromosomes_to_analyze"][wildcards.genome]),
# 		idx = lambda wildcards: expand(
# 			"vcf_joint/{pop}.{gen}.{chrom}.filtered_{ty}.vcf.gz.tbi",
# 			pop=wildcards.population,
# 			gen=wildcards.genome,
# 			ty=wildcards.type,
# 			chrom=config["chromosomes_to_analyze"][wildcards.genome])
# 	output:
# 		"vcf_combined/{population}.{genome}.combined.filtered_{type}.vcf.gz"
# 	params:
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} concat -O z -o {output} {input.vcf}"
#
# rule index_concatenated_vcf:
# 	input:
# 		"vcf_combined/{population}.{gen}.combined.filtered_{type}.vcf.gz"
# 	output:
# 		"vcf_combined/{population}.{gen}.combined.filtered_{type}.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"

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
