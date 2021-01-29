options(stringsAsFactors = FALSE)

library(optparse)


# input: 
# - rpkm_table (transcript)
# - peaks
# - gtf
# - optional parameters
# output:
# - annotated promoters if expressed or not


option_list = list(
  make_option(c("-t", "--rpkm_table"), default=NA, type='character',
              help="RPKM_table. Must contain only desired samples (transcript level. format: tsv)"),
  make_option(c("-p", "--peak_file"), default=NA, type='character',
              help="GRO-seq peaks"),
  make_option(c("-g","--gtf_file"), default=NA, type='character',
              help="Annotation file to identify promoters (default: %default)"),
  make_option(c("-s","--stranded"), default=FALSE, action = 'store_true', type='logical',
              help="Perform TSS-Peak assignment stranded (default: %default)"),
  
  make_option("--distance_cutoff", default=1000, type='integer',
              help="Maximum distance between peak boundaries and TSS (default: %default)"),
  make_option("--expression_threshold", default=3, type='numeric',
              help="Minimum distance between peak boundaries and TSS (default: %default)"),
  make_option("--transcript_size_threshold", default=125, type='integer',
              help="Minimum distance between peak boundaries and TSS (default: %default)"),
  
  make_option("--gene_biotype", default=NULL, type='character',
              help="filter gtf for gene_biotype"),
  
  make_option("--peak_format", default='auto', type='character',
              help="Peak format. Supported \'auto\' \'bed\', \'narrowPeak\'\n\'auto\' is based on suffix (default: %default)"),
  
  make_option(c("-o", "--outprefix"), default=NA, type='character',
              help="output prefix"),
  
  make_option("--export_pairs", default=FALSE, action = 'store_true', type='logical',
              help="Export TSS:Peak pairs")
  
)
opt = parse_args(OptionParser(option_list=option_list))


# Rscript GRO-seq.active_genes.R 
# --export_pairs 
# -r GROseq.featureCounts_tx/GRO-seq.dm6_ensembl96.transcripts.rpkm_scaled.wt_only.tsv 
# -p GROseq.peak_calling/JAMM/peaks.wildtype/peaks/filtered.peaks.narrowPeak 
# -g dm6_ensembl96.gtf 
# -o Stage5_expressed/GROseq_wt.no_dedup.peaks_jamm

# opt$export_pairs= TRUE
# opt$rpkm_table="GROseq.featureCounts_tx/GRO-seq.dm6_ensembl96.transcripts.rpkm_scaled.wt_only.tsv"
# # opt$rpkm_table="GROseq.featureCounts_tx/GRO-seq_umi_dedup.dm6_ensembl96.transcripts.rpkm_scaled.tsv"
# opt$peak_file = "GROseq.peak_calling/JAMM/peaks.wildtype.stranded.narrowPeak"
# # opt$peak_file = "GROseq.peak_calling/MACS2/peaks.wildtype/GRO-seq.wildtype_peaks.narrowPeak"
# opt$gtf_file = "dm6_ensembl96.gtf"
# opt$outprefix = "Stage5_expressed/GROseq_wt.no_dedup.peaks_jamm"
# opt$peak_format = 'auto'
# opt$stranded = TRUE

detect_peak_format = function(peak_file){
  if(grepl('.bed$', peak_file))
    return("bed")
  if(grepl('.narrowPeak$', peak_file))
    return('narrowPeak')
  else
    stop('Not supported')
}


cat('> Loading libraries\n')
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

# choose filtered_bam or bam_umi_dedup here
# rpkm.tab = read_tsv('./featureCounts/GRO-seq.dm6_ensembl96.transcripts.rpkm_scaled.tsv')
# rpkm.tab = read_tsv('./GROseq.featureCounts_tx/GRO-seq_umi_dedup.dm6_ensembl96.transcripts.rpkm_scaled.tsv')
# 
# dm6.tx = subset(import.gff('./dm6_ensembl96.gtf', feature.type = 'transcript'), width > 125)
# f0 = 'GROseq.peak_calling//MACS2/peaks.wildtype/GRO-seq.wildtype_peaks.narrowPeak'
# peaks = import.bed(f0, extraCols = extraCols_narrowPeak)
# 
# f0 = 'GROseq.peak_calling/JAMM/peaks.wildtype/peaks/filtered.peaks.narrowPeak'
# peaks = import.bed(f0, extraCols = extraCols_narrowPeak)
#
# dm6.tss = dm6.tx %>% resize(fix = 'start', width = 1)
# # promoter := TSS + [-300, 500]
# dm6.promoters = dm6.tss %>% resize(fix = 'start', width = 1) %>% resize(fix = 'center', width = 2000)
# start(dm6.promoters) = ifelse(start(dm6.promoters) < 1, 1 ,start(dm6.promoters))
# 
# # has.peak_overlap = dm6.promoters %over% peaks;
# has.peak_overlap = overlapsAny(dm6.promoters, peaks, type = 'any')
cat('> Reading data\n')
dm6.tx = import.gff(opt$gtf_file, feature.type = 'transcript')
if(!is.null(opt$gene_biotype)){
  assertthat::assert_that(all(opt$gene_biotype %in% unique(mcols(dm6.tx)$gene_biotype)), 
                          msg = paste0('Named gene_biotype(s) are not available. Available are', unique(mcols(dm6.tx)$gene_biotype)))
  dm6.tx = subset(dm6.tx, gene_biotype %in% opt$gene_biotype)
}


dm6.tx = subset(dm6.tx, width > opt$transcript_size_threshold)
cat('>> Counting ',length(dm6.tx), 'transcripts\n')

rpkm.tab = read_tsv(opt$rpkm_table)


peaks_file = opt$peak_file
if(opt$peak_format == 'auto'){
  peak_format = detect_peak_format(peaks_file)
} else {
  peak_format = opt$peak_format
}

cat('>> Loading',peak_format,'file format\n')

if(peak_format == 'bed'){
  peaks = import.bed(peaks_file)
} else if(peak_format == 'narrowPeak'){
  peaks = import.bed(peaks_file, extraCols = extraCols_narrowPeak)
} else {
  stop('Peak format not supported')
}

# Identify transcripts TSS with peaks in vicinity (symmetric!)
cat('> Identify TSS with peaks nearby\n')
has.peak = rep(FALSE, length(dm6.tx))
dd0 = distanceToNearest(dm6.tx %>% resize(width = 1, fix='start'), peaks, ignore.strand=!opt$stranded)
has.peak[queryHits(dd0)] = mcols(dd0)$distance < opt$distance_cutoff


# export tss:peak pairs
if(opt$export_pairs){
  cat('>> Exporting pairs\n')
  df.ref = data.frame(dm6.tx, stringsAsFactors = FALSE)[c('seqnames','start','end','strand','transcript_id','gene_id')]
  
  df.add = data.frame(peaks[subjectHits(dd0)] %>% as.data.frame(), 
                      distance = mcols(dd0)$distance,
                      stringsAsFactors = FALSE)
  
  dfout = cbind(df.ref[queryHits(dd0),], df.add)
  i.not = setdiff(1:length(dm6.tx), queryHits(dd0))

  dummy.tab = matrix(NA, nrow = length(i.not), ncol = ncol(dfout) - ncol(df.ref))
  colnames(dummy.tab) = colnames(dfout)[(ncol(df.ref)+1):ncol(dfout)]
  dfout = rbind(dfout, 
                cbind(df.ref[i.not,], dummy.tab))
  print(nrow(dfout))
  # print(head(dfout))
  
  write_tsv(dfout, paste0(opt$outprefix,'.closest_peak_mapping.tsv'))
}

# Identify transcripts that are sufficiently read
cat('> Identify expressed transcripts\n')
rpkm.tab_rowMeans = rowMeans(rpkm.tab[,-1])
txid2rpkm_mean = cbind(rpkm.tab[,1],rpkm.tab_rowMeans) %>% deframe
is.expressed = txid2rpkm_mean[mcols(dm6.tx)$transcript_id] > opt$expression_threshold

# Annotated transcripts
cat('> Annotating transcripts\n')
mcols(dm6.tx)$peak_distance = NA
mcols(dm6.tx)$peak_distance[queryHits(dd0)] = mcols(dd0)$distance
mcols(dm6.tx)$peak_qValue = NA
mcols(dm6.tx)$peak_qValue[queryHits(dd0)] = mcols(peaks)$qValue[subjectHits(dd0)]
mcols(dm6.tx)$peak_signalValue = NA
mcols(dm6.tx)$peak_signalValue[queryHits(dd0)] = mcols(peaks)$signalValue[subjectHits(dd0)]
mcols(dm6.tx)$peak_strand = NA
mcols(dm6.tx)$peak_strand[queryHits(dd0)] = strand(peaks)[subjectHits(dd0)]

mcols(dm6.tx)$expression = txid2rpkm_mean[mcols(dm6.tx)$transcript_id]

mcols(dm6.tx)$target_expressed = is.expressed
mcols(dm6.tx)$target_peak = has.peak
# mcols(dm6.tx)$target = is.expressed & has.peak

score(dm6.tx) = txid2rpkm_mean[mcols(dm6.tx)$transcript_id]
names(dm6.tx) = mcols(dm6.tx)$transcript_id

df.stats = data.frame(
  transcript_id = mcols(dm6.tx)$transcript_id, 
  gene_id = mcols(dm6.tx)$gene_id,
  peaks = has.peak,
  expressed = is.expressed,
  is.active = has.peak & is.expressed)

df.promoters = df.stats %>% 
  summarize(n.peaks = sum(has.peak), 
            n.expressed = sum(is.expressed), 
            n.active = sum(has.peak&is.expressed))

df.genes = df.stats %>% select(-transcript_id) %>% 
  group_by(gene_id) %>% 
  summarize(peaks = any(peaks, na.rm = TRUE),
            expressed = any(expressed, na.rm = TRUE),
            is.active = any(is.active, na.rm = TRUE)) %>% 
  summarize(n.peaks = sum(peaks), 
            n.expressed = sum(expressed), 
            n.active = sum(is.active))

df.stats = rbind(cbind('type' = 'transcripts' ,df.promoters, total = nrow(df.stats)), 
                 cbind('type' = 'genes', df.genes, total = length(df.stats$gene_id %>% unique())))

cat('> Writing tables\n')
write_tsv(df.stats, paste0(opt$outprefix, '.stats.tsv'))


dm6.tx %>% as.data.frame() %>% head
dm6.anno = dm6.tx %>% as.data.frame() %>% 
  dplyr::select('seqnames','start','end','strand','transcript_id','gene_id','gene_name',
                'peak_distance','peak_qValue','peak_signalValue','peak_strand','expression',
                'target_expressed','target_peak') #, 'target')
write_tsv(dm6.anno, paste0(opt$outprefix,'.annotated.tsv'))

gr0 = subset(dm6.tx, target_peak)
export.bed(gr0, paste0(opt$outprefix, '.peak.transcripts.bed'))
export.bed(gr0 %>% resize(fix = 'start', width =1) %>% unique(),
           paste0(opt$outprefix, '.peak.TSS_uniq.bed'))

gr0 = subset(dm6.tx, !target_peak)
export.bed(gr0, paste0(opt$outprefix, '.no_peak.transcripts.bed'))
export.bed(gr0 %>% resize(fix = 'start', width =1) %>% unique(),
           paste0(opt$outprefix, '.no_peak.TSS_uniq.bed'))


cat('> ...Done...\n')
