
rule featureCounts:
    input:
        bams=expand(os.path.join("filtered_bam","{sample}.bam"), sample = samplenames),
        annotation=config['annotation_gtf']
    output:
        table=os.path.join("featureCounts","GRO-seq.featureCounts.mapq{mapq}.tsv")
    params:
        # mapq filter?
        # multimappers?
        feature="-t gene",
        multimappers="", # not --fraction
        mapq_filter="-Q {mapq}",
        stranded="-s 1"
    log: os.path.join("featureCounts","log","GRO-seq.featureCounts.mapq{mapq}.log")
    threads: 16
    shell:
        "featureCounts -T {threads} {params} -a {input.annotation} -o {output.table} {input.bams} &> {log}"


rule featureCounts_dedup:
    input:
        bams=expand(os.path.join("bam_umi_dedup","{sample}.bam"), sample = samplenames),
        annotation=config['annotation_gtf']
    output:
        table=os.path.join("featureCounts_dedup","GRO-seq.featureCounts.mapq{mapq}.tsv")
    params:
        # mapq filter?
        # multimappers?
        feature="-t gene",
        multimappers="", # not --fraction
        mapq_filter="-Q {mapq}",
        stranded="-s 1"
    log: os.path.join("featureCounts_dedup","log","GRO-seq.featureCounts.mapq{mapq}.log")
    threads: 16
    shell:
        "featureCounts -T {threads} {params} -a {input.annotation} -o {output.table} {input.bams} &> {log}"
