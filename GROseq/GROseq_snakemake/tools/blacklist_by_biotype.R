suppressPackageStartupMessages(library(optparse))


# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list(
  make_option(c("-g", "--reference"), type="character", 
              help="Genome reference (gtf)"),
  make_option(c('-i',"--fasta_index"),  type = 'character',
              help = "Fasta index to check chromosome boundaries"),
  
  make_option("--flank", default=50, type = 'integer',
              help=" Flank to add to blacklisted genes [default: %default]"),
  make_option("--biotypes", default="rRNA,snRNA,snoRNA,tRNA", 
              help=" Biotypes to black list. Comma-separated [default: %default]"),

  make_option(c('-o',"--outfile"), type = 'character',
              help="Write blacklisted genes to this file (bed format)")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

# opt = list()
# opt$gtf = './dm3_ensembl78.gtf'
# opt$faindex <- './dm3.fa.idx'
# opt$biotypes = c('rRNA','snRNA','snoRNA','tRNA')
# opt$flank = 100


suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

chr.converter = deframe(cbind('dmel_mitochondrion_genome', 'mitochondrion_genome'))

convert_chromosomes <- function(chrom.set, converter){
  b <- chrom.set %in% names(converter)
  if(length(b) > 0){
    cat('Converting chromosomes:', chrom.set[b],'\n')
    chrom.set[b] = converter[chrom.set[b]]
  }
  return(chrom.set)
}

ref.genes <- import.gff(opt$reference, feature.type  ='gene')
chrom.size <- deframe(read_tsv(opt$fasta_index, col_names = FALSE)[,1:2])
biotypes = strsplit(opt$biotypes,',')[[1]]

seqlevels(ref.genes) <- convert_chromosomes(seqlevels(ref.genes), chr.converter)
names(chrom.size) <- convert_chromosomes(names(chrom.size), chr.converter)

assertthat::assert_that(all(seqlevels(ref.genes) %in% names(chrom.size)), 
                        msg = paste('Missing chromomsomes in index', seqlevels(ref.genes)[!seqlevels(ref.genes) %in% names(chrom.size)]))

cat("filtering biotypes:", opt$biotypes, '\n')

ref.blacklist <- subset(ref.genes, gene_biotype %in% biotypes)
ref.blacklist_ext <- resize(ref.blacklist, width = width(ref.blacklist) + opt$flank * 2, fix = 'center')


ii <- which(start(ref.blacklist_ext) < 0 | end(ref.blacklist_ext) > chrom.size[as.vector(seqnames(ref.blacklist_ext))])
cat("Counting genes that exceed chromosome boundaries:", length(ii),'\n')


if(length(ii) > 0){
  cat('Fixing to match boundaries...\n')

  start(ref.blacklist_ext)[start(ref.blacklist_ext) < 1] = 1
  i <- end(ref.blacklist_ext) > chrom.size[as.vector(seqnames(ref.blacklist_ext))]
  end(ref.blacklist_ext)[i] = chrom.size[as.vector(seqnames(ref.blacklist_ext))][i]
} else {
  cat('No flanks to fix\n')
}

score(ref.blacklist_ext) <- 0
names(ref.blacklist_ext) <- mcols(ref.blacklist_ext)$gene_id

export.bed(object = ref.blacklist_ext, con = opt$outfile)

