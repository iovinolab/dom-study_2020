# require bedtool2
localrules: link_annotationgtf, extract_genes

rule link_annotationgtf:
    input:
        gtf=config['annotation_gtf']
    output:
        gtf=os.path.join("Annotations","annotation.gtf")
    shell:
        "cp {input.gtf} {output.gtf}"

rule extract_genes:
    input:
        gtf=rules.link_annotationgtf.output.gtf,
        chroms_ignore = config['chroms_ignore']
    output:
        bed=os.path.join("Annotations","genes.bed")
    params:
        chroms_ignore = config['chroms_ignore']
    shell:
        """
        cat {input.gtf} | awk '$3 == "gene"' | egrep -v $(paste -sd '|' {input.chroms_ignore}) | gff2bed > {output.bed}
        """

rule extract_transcript:
    input:
        gtf=rules.link_annotationgtf.output.gtf,
        chroms_ignore = config['chroms_ignore']
    output:
        transcripts=os.path.join("Annotations","transcripts.bed")
    params:
        chroms_ignore = config['chroms_ignore']
    shell:
        """
        cat {input.gtf} | awk '$3 == "transcript"' | egrep -v $(paste -sd '|' {input.chroms_ignore}) | gff2bed > {output.transcripts}
        """

rule extract_tss:
    input:
        gtf=rules.extract_transcript.output.transcripts,
        genome_index=config['genome_index']
    output:
        bed=os.path.join("Annotations","transcript_start_sites.bed")
    params:
        params="-s -l 1"
    shell:
        """
        bedtools flank -s -l 50 -r 50 -i {input.gtf} -g {input.genome_index} > {output.bed}
        """
