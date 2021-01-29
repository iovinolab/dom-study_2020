library(optparse)

chroms.valid = paste0(c('2L','2R','3L','3R','4','X'), collapse =',')

option_list <- list( 
  make_option(c("-g", "--gtf"), 
              help="Genome annotation file (gtf)"),
  make_option(c("-o", "--outfile_yaml"), 
              help="Metafile for output (yaml)"),
  
  make_option("--select_chromosomes", default = chroms.valid,
              help="Limit segmentation to given chromosomes (default: %default)"),
  
  make_option(c("--lib.location"), default = NULL, type = 'character',
              help="Specify library location")
  
)
opt <- parse_args(OptionParser(option_list=option_list))

lib.loc = opt$lib.location

chroms.valid <- strsplit(opt$select_chromosomes,',')[[1]]

# suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer, lib.loc = lib.loc))
suppressPackageStartupMessages(library(dplyr, lib.loc = lib.loc))
suppressPackageStartupMessages(library(tibble, lib.loc = lib.loc))
suppressPackageStartupMessages(library(yaml, lib.loc = lib.loc))



# txdb0 <- makeTxDbFromGFF('/data/repository/organisms/dm3_ensembl/ensembl/release-78/genes.gtf')
anno.transcript <- import.gff(opt$gtf, feature.type = 'transcript')
anno.transcript <- subset(anno.transcript, seqnames %in% chroms.valid)
tx2gid <- deframe(mcols(anno.transcript)[c('transcript_id','gene_id')])

tx_set <- anno.transcript %>% unique %>% as.data.frame()
tss.pos <- tx_set %>% 
  filter(strand == '+') %>% 
  mutate(start = start,
         end = start + 200) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

tss.neg <- tx_set %>% 
  filter(strand == '-') %>% 
  mutate(start = end - 200,
         end = end) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


tss_gr <- sort(c(tss.pos, tss.neg))

tx_unit.pos <- tx_set %>% 
  filter(strand == '+') %>% 
  mutate(start = ifelse(width >= 600, start + 400, start),
         end = end) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

tx_unit.neg <- tx_set %>% 
  filter(strand == '-') %>% 
  mutate(start = start,
         end = ifelse(width >= 600, end - 400, end)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

tx_unit_gr <- sort(c(tx_unit.pos, tx_unit.neg))


tss_ds.pos <- tx_set %>% 
  filter(strand == '+') %>% 
  mutate(start = start + 201,
         end = start + 401) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

tss_ds.neg <- tx_set %>% 
  filter(strand == '-') %>% 
  mutate(start = end - 401,
         end = end - 201) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

tss_ds_gr <- sort(c(tss_ds.pos, tss_ds.neg))


prefix = gsub('(.*)\\.(yaml|yml)$','\\1',opt$outfile_yaml)
# prefix=paste0(opt$outdir, gsub('(.*)\\.(gtf|gff)','\\1',basename(opt$gtf)))

df0 <- data.frame(tss = paste0(prefix,'.segment_tss.gff'),
  tssds = paste0(prefix,'.segment_tssds.gff'),
  txunit = paste0(prefix, '.segment_txunit.gff'),
  gene_body = opt$gtf, stringsAsFactors = FALSE)

yaml::write_yaml(df0, opt$outfile_yaml)

export.gff3(tss_gr, df0$tss)
export.gff3(tx_unit_gr, df0$txunit)
export.gff3(tss_ds_gr, df0$tssds)
