rule computeMatrix_tss:
    input:
        bigwigs=expand(os.path.join("bamCoverage", "{sample}.{bam_folder}.CPM_mapq.{strand}.bw"),
            bam_folder = "{dedup}",sample = samplenames, strand = 'both'),
        gtf=config['annotation_gtf']
    output:
        matrix=os.path.join("deeptools_visualization", "computeMatrix_genes.{dedup}.mat.gz")
    params:
        referencepoint="--referencePoint TSS",
        window="-b 250 -a 1000",
        binsize="--binSize 25",
        missingdata="--missingDataAsZero"
    threads: 24
    params:
        labals="--samplesLabels {names}".format(names = samplenames)
    shell:
        "computeMatrix reference-point -p {threads} {params} -S {input.bigwigs} -R {input.gtf} -o {output.matrix}"


rule heatmap_tss:
    input:
        matrix=rules.computeMatrix_tss.output.matrix
    output:
        png=os.path.join("deeptools_visualization", "heatmap_genes.{dedup}.png"),
        regions=os.path.join("deeptools_visualization", "heatmap_genes.{dedup}.sorted.bed"),
    shell:
        "plotHeatmap -m {input.matrix} -o {output.png} --outFileSortedRegions {output.regions}"
