library(optparse)

option_list <- list(
  make_option(c("-g", "--gtf"), 
              help="Genome annotation file (gtf)"),
  make_option(c("-p", "--promoter"),
              help="Quantification in promoter region"),
  make_option(c("-d", "--downstream_feature"),
              help="Quantification in promoter region"),

  make_option(c('-o','--outprefix'),
              help = 'output prefix for GROseq pausing tables per sample (tsv)'),
  
  make_option(c("--lib.location"), default = NULL, type = 'character',
              help="Specify library location")
)


opt <- parse_args(OptionParser(option_list=option_list))

lib.loc = opt$lib.location

suppressPackageStartupMessages(library(rtracklayer, lib.loc = lib.loc))
suppressPackageStartupMessages(library(readr, lib.loc = lib.loc))
suppressPackageStartupMessages(library(dplyr, lib.loc = lib.loc))
suppressPackageStartupMessages(library(tibble, lib.loc = lib.loc))
suppressPackageStartupMessages(library(reshape2, lib.loc = lib.loc))
suppressPackageStartupMessages(library(ggplot2, lib.loc = lib.loc))

## debug
# opt$gtf = 'dm6_ensembl96.gtf'
# opt$promoter = 'pausing.bam_umi_dedup/test_run.tss.counts.tsv'
# opt$downstream_feature = 'pausing.bam_umi_dedup/test_run.tssds.counts.tsv'
#
# Rscript GROseq_pausing.R
#  -g dm6_ensembl96.gtf 
#  -p pausing.bam_umi_dedup/test_run.tss.counts.tsv 
#  -d pausing.bam_umi_dedup/test_run.tssds.counts.tsv
#  -o pausing.bam_umi_dedup/test_run


anno.transcript = import.gff(opt$gtf, feature.type = 'transcript')
tx2gid <- deframe(mcols(anno.transcript)[c('transcript_id','gene_id')])
tx2biotype = deframe(mcols(anno.transcript)[c('transcript_id','transcript_biotype')])

anno.genes <- import.gff(opt$gtf, feature.type = 'gene')
gid2biotype = deframe(mcols(anno.genes)[c('gene_id','gene_biotype')])
gid2symbol = deframe(mcols(anno.genes)[c('gene_id','gene_name')])
gid2width <- deframe(data.frame(mcols(anno.genes)$gene_id, width(anno.genes), stringsAsFactors = FALSE))

tabset <- list(tss = opt$promoter, txunit = opt$downstream_feature)

quantTabSet <- lapply(tabset, read_tsv, col_names = TRUE, comment = '#')
# quantTabSet <- lapply(quantTabSet, function(x) {colnames(x) <- gsub('.*(T.*)\\.bam$','\\1',colnames(x));return(x)})
quantTabSet <- lapply(quantTabSet, function(x) {colnames(x) <- gsub('.bam$','',basename(colnames(x)));return(x)})


cols.anno = c('Geneid','Chr','Start','End','Strand','Length')
dd0 <- melt(quantTabSet, id.vars = cols.anno)
colnames(dd0)[1] = 'Transcriptid'

head(dd0)

dd.wide <- dcast(dd0, formula(paste(paste(c('Transcriptid','variable'), collapse = '+')," ~ L1")))
# Active definition (active/inactive):
# - GRO-seq in TSS or TXunit
# Active TSS (active/inactive)
# - GRO-seq in TSS 
# Paused (gradual, 3 groups) [Inf - 1, 1-0.5, 0.5 - -Inf]
# - log2 TSS/TXUNIT
# unique(sort(dd.wide$txunit)) %>% head
# ggplot(dd.wide, aes(tss)) + geom_density() + geom_vline(xintercept = 1e-2) + scale_x_log10() + facet_wrap(~condition)
# ggplot(dd.wide, aes(txunit)) + geom_density() + scale_x_log10() + facet_wrap(~condition)

dd.wide <- dd.wide %>% 
  mutate(pausing.ratio = (tss + 1e-8)/(txunit + 1e-8)) %>% 
  mutate(Geneid = tx2gid[Transcriptid],
 # Definition of 'active' run by RPKM + peak calling. Not necessary any more in this part 
         tss.active = tss > 1e-2) %>% 
  ## Be stricter with transcript_active
  ## reconsider txunit > 1e-3. Upstream TSS may be considered active for downstream TSS beig txunit
  mutate(transcript_active = tss > 1e-2 | txunit > 1e-2) %>% 
  # group by signal in TSS
  group_by(Geneid) %>% 
  mutate(gene_active = any(transcript_active)) %>% 
  # not.paused := log2(pausing) <= 0
  mutate(pausing.status = ifelse(transcript_active & log2(pausing.ratio) > 0.5, 'paused',   
                         # relaxed cutoff - e.g. 0.5    
                         ifelse(transcript_active & log2(pausing.ratio) > 1,'strongly.paused', 'not.paused')), 
         symbol = gid2symbol[Geneid]) %>% 
  ungroup

dd.wide_per_sample = split(dd.wide %>% select(-variable), dd.wide$variable)

for(sample0 in names(dd.wide_per_sample)){
  write_tsv(dd.wide_per_sample[[sample0]], paste0(opt$outprefix,'.',sample0,'.tsv'))
}
