# requires deeptools

# First strand synthesis only, thus reverse for deeptools application
strand_map = {'forward': 'reverse','reverse': 'forward'}

rule deeptools_plotFingerprint:
    input:
        bams=expand(os.path.join("{dedup}", "{sample}.bam"),
            dedup = "{dedup}",sample = samplenames),
        bais=expand(os.path.join("{dedup}", "{sample}.bam.bai"),
            dedup = "{dedup}",sample = samplenames),
    output:
        png=os.path.join("deeptools","plotFingerprint.{dedup}.png"),
        counts=os.path.join("deeptools","plotFingerprint.{dedup}.counts.tsv"),
        metrics=os.path.join("deeptools","plotFingerprint.{dedup}.qc_metric.tsv")
    params:
        labels=" ".join(samplenames)
    log:
        os.path.join("deeptools","log","plotFingerprint.{dedup}.log")
    threads: 16
    shell:
        """
        plotFingerprint -p {threads} --labels {params.labels} -b {input.bams} -o {output.png} \
                --outRawCounts {output.counts} \
                --outQualityMetrics {output.metrics} \
                &> {log}
        """

rule deeptools_plotEnrichment:
    input:
        bams=expand(os.path.join("{dedup}", "{sample}.bam"),
            dedup = "{dedup}",sample = samplenames),
        bais=expand(os.path.join("{dedup}", "{sample}.bam.bai"),
            dedup = "{dedup}",sample = samplenames),
        gtf=config['annotation_gtf']
    output:
        png=os.path.join("deeptools","plotEnrichment.{dedup}.png")
    log:
        os.path.join("deeptools","log","plotEnrichment.{dedup}.log")
    params:
        persample="--perSample"
    threads:18
    shell:
        """
        plotEnrichment -p {threads} --bamfiles {input.bams} --BED {input.gtf} --plotFile {output.png} &> {log}
        """

rule deeptools_estimateReadFilter:
    input:
        bams=expand(os.path.join("{dedup}", "{sample}.bam"),
            dedup = "{dedup}",sample = samplenames),
        bais=expand(os.path.join("{dedup}", "{sample}.bam.bai"),
            dedup = "{dedup}",sample = samplenames)
    output:
        table=os.path.join("deeptools","estimateReadFiltering.{dedup}.tsv")
    params:
        mapq_filter="--minMappingQuality {mapq}".format(mapq = config['mapq_filter'])
    threads: 8
    shell:
        "estimateReadFiltering -p {threads} {params} -b {input.bams} -o {output.table}"


rule deeptools_multiBigwigSummary:
    input:
        bigwigs=expand(os.path.join("bamCoverage","{sample}.{bam_folder}.CPM_mapq.{strand}.bw"),
            sample = samplenames, bam_folder = 'bam_umi_dedup', strand = 'both'),
        bed=rules.extract_genes.output.bed
    output:
        npz=os.path.join("deeptools","multiBigwigSummary.genes.{dedup}.npz"),
        counts=os.path.join("deeptools","multiBigwigSummary.genes.{dedup}.counts.tsv")
    log:
        os.path.join("deeptools","log","multiBigwigSummary.genes.{dedup}.log"),
    threads: 16
    shell:
        """
        multiBigwigSummary BED-file -p {threads} -b {input.bigwigs} --BED {input.bed} -o {output.npz} --outRawCounts {output.counts} &> {log}
        """

rule deeptools_plotCorrelation:
    input:
        npz=rules.deeptools_multiBigwigSummary.output.npz
    output:
        png=os.path.join("deeptools","plotCorrelation.{dedup}.png"),
        matrix=os.path.join("deeptools","plotCorrelation.{dedup}.correlations.tsv"),
    params:
        plot="--whatToPlot heatmap",
        corMethod="--corMethod pearson",
        others="--skipZeros --removeOutliers --plotNumbers",
        mapq="--minMappingQuality {mapq}".format(mapq = config['mapq_filter'])
    log:
        os.path.join("deeptools","log","plotCorrelation.{dedup}.log")
    shell:
        """
        plotCorrelation {params.plot} {params.corMethod} {params.others} --corData {input.npz} --outFileCorMatrix {output.matrix} -o {output.png} &> {log}
        """

rule deeptools_plotPCA:
    input:
        npz=rules.deeptools_multiBigwigSummary.output.npz
    output:
        png=os.path.join("deeptools","plotPCA.{dedup}.png"),
        matrix=os.path.join("deeptools","plotPCA.{dedup}.counts.tsv")
    log:
        os.path.join("deeptools","log","plotPCA.{dedup}.log")
    shell:
        """
        plotPCA --transpose --corData {input.npz} --outFileNameData {output.matrix} -o {output.png} &> {log}
        """

rule deeptools_bamCoverage:
    input:
        bam=os.path.join("{bam_folder}","{sample}.bam"),
        bai=os.path.join("{bam_folder}","{sample}.bam.bai")
    output:
        bigwig=os.path.join("bamCoverage","{sample}.{bam_folder}.coverage_raw.{strand}.bw")
    params:
        filterRNAstrand=lambda wildcards: "" if not wildcards.strand in strand_map.keys() else ' '.join(["--filterRNAstrand ", strand_map[wildcards.strand]]),
        # Offset="--Offset 1",
        binsize="--binSize 1"
    threads: 8
    log:
        os.path.join("bamCoverage","log","{sample}.{bam_folder}.coverage_raw.{strand}.log")
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} {params} &> {log}"

ruleorder: deeptools_bamCoverage_mapq > deeptools_bamCoverage_norm

rule deeptools_bamCoverage_mapq:
    input:
        bam=os.path.join("{bam_folder}","{sample}.bam"),
        bai=os.path.join("{bam_folder}","{sample}.bam.bai")
    output:
        bigwig=os.path.join("bamCoverage","{sample}.{bam_folder}.coverage_mapq.{strand}.bw")
    params:
        filterRNAstrand=lambda wildcards: "" if not wildcards.strand in strand_map.keys() else ' '.join(["--filterRNAstrand ", strand_map[wildcards.strand]]),
        # Offset="--Offset 1",
        binsize="--binSize 1",
        mapq="--minMappingQuality {mapq}".format(mapq = config['mapq_filter']),

    threads: 8
    log:
        os.path.join("bamCoverage","log","{sample}.{bam_folder}.coverage_mapq.{strand}.log")
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} {params} &> {log}"


rule deeptools_bamCoverage_norm:
    input:
        bam=os.path.join("{bam_folder}","{sample}.bam"),
        bai=os.path.join("{bam_folder}","{sample}.bam.bai")
    output:
        bigwig=os.path.join("bamCoverage","{sample}.{bam_folder}.{normalization}_mapq.{strand}.bw")
    params:
        filterRNAstrand=lambda wildcards: "" if not wildcards.strand in strand_map.keys() else ' '.join(["--filterRNAstrand ", strand_map[wildcards.strand]]),
        normalize=lambda wildcards: "--normalizeUsing {norm} --effectiveGenomeSize {genome_size}".format(norm = wildcards.normalization, genome_size = str(config['genome_size'])) if wildcards.normalization == "RPGC" else "--normalizeUsing {norm}".format(norm = wildcards.normalization),
        binsize="--binSize 1",
        mapq="--minMappingQuality {mapq}".format(mapq = config['mapq_filter']),
        chroms_ignore="--ignoreForNormalization $(cat {chrom_file} | paste -sd ' ')".format(chrom_file = config['chroms_ignore'])
    threads: 8
    log:
        os.path.join("bamCoverage","log","{sample}.{bam_folder}.{normalization}_mapq.{strand}.log")
    shell:
        "bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} {params} &> {log}"
