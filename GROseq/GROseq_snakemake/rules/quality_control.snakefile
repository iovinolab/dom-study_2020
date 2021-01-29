
# indirs = multiqc_input_check(return_value = "indir")
# multiqc_input_check(return_value = "infiles")
localrules: multiqc, multiqc_fastq

rule fastqc:
    input:
        fastq=os.path.join("{folder}","{sample}.fastq.gz")
    output:
        fastqc_rep=os.path.join("{folder}","{sample}_fastqc.html")
    threads: 4
    shell:
        "fastqc -t {threads} {input.fastq}"

rule multiqc:
    input:
        expand(os.path.join("fastq_trimmed", "{sample}.fastq.gz"), sample = samplenames),
            expand(os.path.join("Bowtie2", "{sample}.bam"), sample = samplenames),
            expand(os.path.join("bam_umi_dedup", "{sample}.bam"), sample = samplenames),
            # rules.deeptools_multiBamSummary.output.counts,
            # rules.deeptools_plotEnrichment.output.png,
            # rules.deeptools_plotFingerprint.output.png,
            # rules.deeptools_estimateReadFilter.output.table,
            # rules.deeptools_plotCorrelation.output.matrix,
            # rules.deeptools_plotPCA.output.png
            expand(os.path.join("{bam_folder}","{sample}.idxstats"), bam_folder = 'bam_umi_dedup', sample = samplenames),
            expand(os.path.join("deeptools", "estimateReadFiltering.{dedup}.tsv"), dedup = 'Bowtie2'),
            expand(os.path.join("deeptools", "multiBigwigSummary.genes.{dedup}.counts.tsv"), dedup = "bam_umi_dedup"),
            expand(os.path.join("deeptools", "plotCorrelation.{dedup}.correlations.tsv"), dedup = "bam_umi_dedup"),
            expand(os.path.join("deeptools", "plotFingerprint.{dedup}.counts.tsv"), dedup = "bam_umi_dedup"),
            expand(os.path.join("deeptools", "plotFingerprint.{dedup}.qc_metric.tsv"), dedup = "bam_umi_dedup")
    output:
        os.path.join("Multiqc","multiqc_report.html")
    params:
        directories=' '.join(["fastq_trimmed", "Bowtie2", "bam_umi_dedup", "deeptools"]),
        directory_out=directory(os.path.join("Multiqc"))
    log:
        os.path.join("Multiqc","log","multiqc.log")
    shell:
        "multiqc -s -f --outdir {params.directory_out} {params.directories} &> {log}"

rule multiqc_fastq:
    input:
        files=expand(os.path.join("Fastq","{sample}_fastqc.html"), sample = samplenames)
    output:
        multiqc=os.path.join('Fastq_multiqc','multiqc_report.html')
    params:
        directory="Fastq",
        directory_out="Fastq_multiqc"
    log:
        os.path.join("Fastq_multiqc","log","multiqc.log")
    shell:
        "multiqc -f -s --outdir {params.directory_out} {params.directory} &> {log}"



##
## snakePipes adoption
##
# rule multiQC:
#     input:
#         multiqc_input_check(return_value = "infiles")
#     output: "multiQC/multiqc_report.html"
#     params:
#         indirs = multiqc_input_check(return_value = "indir")
#     log:
#         out = "multiQC/multiQC.out",
#         err = "multiQC/multiQC.err"
#     conda: CONDA_SHARED_ENV
#     shell:
#         "multiqc -o multiQC -f {params.indirs} > {log.out} 2> {log.err}"
