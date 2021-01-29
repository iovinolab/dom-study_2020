rule pausing:
    input:
        annotation_gtf=config['annotation_gtf'],
        samplesheet=config['samplesheet'],
        bams=expand(os.path.join('{section}','filtered_bam','{sample}.bam'), sample = samplenames, section = "{section}")
    output:
        done=os.path.join("{section}",'Pausing','GROseq_pausing.DONE')
    params:
        script=os.path.join(maindir,'tools','GROseq_pausing.sh'),
        prefix=os.path.join("{section}",'Pausing','GROseq_pausing'),
    log: os.path.join("{section}",'Pausing','GROseq_pausing','log','GROseq_pausig.log')
    threads: 12
    shell:
        "bash {params.script} {threads} {params.prefix} {input.annotation_gtf} {input.samplesheet} {input.bams} &> {log}"
