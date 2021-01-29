# requires sambamba bowtie2 umi_tools seqtk FastQC

# ruleorder: Bowtie2 > filter_bam > index > flagstat # > dedup
#
ruleorder: filter_bam > index
localrules: link_fastq, copy_config, idxstats

# requires  cutadapt umi_tools subread

rule link_fastq:
    input:
        fastq=os.path.join(config["fastq_directory"],"{sample}.fastq.gz")
    output:
        fastq=os.path.join("Fastq","{sample}.fastq.gz")
    shell:
        "ln -fs {input.fastq} {output.fastq}"

rule copy_config:
    output:
        yaml="config.yaml"
    run:
        yaml_stream = open(output.yaml,"w")
        yaml.dump(config, yaml_stream)
        yaml_stream.close()

rule trim_reads:
    input:
        fastq=os.path.join("Fastq","{sample}.fastq.gz")
    output:
        fastq=temp(os.path.join("fastq_trimmed","{sample}_trimmed.fastq.gz"))
    params:
        adapter_seq="N\{4\}" + config["adapter_sequence"],
        # need to replace --cut. It's applied before adapter trimming...
        params="--overlap=3 --minimum-length=12 --max-n 0" # default: --overlap 3
    log:
        os.path.join("fastq_trimmed","log","{sample}_cutadapt.log")
    threads: 8
    shell:
        """
        cutadapt -a {params.adapter_seq} -j {threads} {params.params} -o {output.fastq} {input.fastq} &> {log}
        """

rule umi_extract:
    input:
        fastq=rules.trim_reads.output.fastq
    output:
        fastq=os.path.join("fastq_trimmed","{sample}.fastq.gz")
    params:
        bcpattern="NNNN"
    log:
        os.path.join("fastq_trimmed","log","{sample}_umi_extract.log")
    shell:
        "umi_tools extract --bc-pattern={params.bcpattern} -I {input.fastq} -S {output.fastq} -L {log}"


rule Bowtie2:
    input:
        fastq=os.path.join("fastq_trimmed","{sample}.fastq.gz"),
        index=os.path.join(config["bowtie2_index"],"genome.fa")
    output:
        bam=os.path.join("Bowtie2", "{sample}.bam"),
        bam_tmp=temp(os.path.join("Bowtie2","{sample}_bowtie2.bam")),
        metrics=os.path.join("Bowtie2","metrics" ,"{sample}.metrics.tsv")
    params:
        index=os.path.join(config["bowtie2_index"],"genome"),
        alignment_mode="--local"
    threads: 8
    log:
        error=os.path.join("Bowtie2", "log", "{sample}.errlog"),
    shell:
        """
        bowtie2 -p {threads} --met-file {output.metrics} {params.alignment_mode} -x {params.index} -U {input.fastq} 2> {log.error} |
            sambamba -q view -t {threads} -S -f bam -o {output.bam_tmp} /dev/stdin &&
            sambamba -q sort -t {threads} -o {output.bam} {output.bam_tmp} && sambamba index {output.bam}
        """

# should collect all read level filterings
## - blacklists
## - unwanted chromosomes (rDNA?)
## - No MAPQ, that's for downstream programs
rule filter_bam:
    input:
        bam=rules.Bowtie2.output.bam
    output:
        bam=os.path.join("filtered_bam","{sample}.bam"),
        bai=os.path.join("filtered_bam","{sample}.bam.bai")
    params:
        mapq=config['mapq_filter']
    threads: 8
    log:
        os.path.join("filtered_bam","log","{sample}.log"),
    shell:
        "sambamba -q view -f bam -t {threads} -o {output.bam} -F \"mapping_quality >= {params.mapq}\" {input.bam}"

rule umi_dedup:
    input:
        bam=rules.filter_bam.output.bam
    output:
        bam=os.path.join("bam_umi_dedup","{sample}.bam"),
    params:
        stats=os.path.join("bam_umi_dedup","{sample}")
    log:
        os.path.join("bam_umi_dedup","log","{sample}.log")
    shell:
        "umi_tools dedup --output-stats {params.stats} -I {input.bam} -S {output.bam} &> {log}"

rule index:
    input:
        bam=os.path.join("{bam_folder}","{sample}.bam")
    output:
        index=os.path.join("{bam_folder}","{sample}.bam.bai")
    log:
        os.path.join("{bam_folder}","log","{sample}_index.log")
    threads: 8
    shell:
        "sambamba -q index -t {threads} {input.bam} {output.index}"

rule flagstat:
    input:
        bam=os.path.join("{bam_folder}","{sample}.bam")
    output:
        flagstat=os.path.join("{bam_folder}","{sample}.flagstat")
    log:
        os.path.join("{bam_folder}","log","{sample}_flagstat.log")
    threads: 8
    shell:
        "sambamba -q flagstat -t {threads} {input.bam} 1> {output.flagstat} 2> {log}"

rule idxstats:
    input:
        bam=os.path.join("{bam_folder}","{sample}.bam")
    output:
        idxstats=os.path.join("{bam_folder}","{sample}.idxstats")
    log:
        os.path.join("{bam_folder}","log","{sample}_idxstats.log")
    shell:
        "samtools idxstats {input.bam} 1> {output.idxstats} 2> {log}"
