library(optparse)

option_list <- list( 
  make_option(c("-c", "--counts_feature"), 
              help="featureCounts of gene feature"),
  make_option(c("-t", "--counts_total"), 
              help="featureCounts of gene"),
  make_option(c("-o", "--outfile"), 
              help="Output of normalization"),
  
  make_option(c("--lib.location"), default = NULL, type = 'character',
              help="Specify library location")
)
opt <- parse_args(OptionParser(option_list=option_list))

lib.loc = opt$lib.location

library(readr, lib.loc = lib.loc)
library(dplyr, lib.loc = lib.loc)
library(tibble, lib.loc = lib.loc)
library(reshape2, lib.loc = lib.loc)
library(ggplot2, lib.loc = lib.loc)

# michael's version
# https://support.bioconductor.org/p/91218/
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

tpm4 <- function(counts, libsize, len){
  x <- counts/len
  return(t(t(x)*1e6/libsize))
}

genes.tab <- read_tsv(opt$counts_total, comment=  '#')
cts.tab <- read_tsv(opt$counts_feature, comment = '#')

libsize0 <- colSums(genes.tab[colnames(cts.tab)[-c(1:6)]])
cts.norm <- cbind(cts.tab[,1:6], tpm4(cts.tab[-c(1:6)], libsize = libsize0 , len = cts.tab$Length))
write_tsv(x = cts.norm, opt$outfile)


# genes.tab <- read_tsv('../data/Pol2pausing/create_pausing_regions.test/featureCounts_genes.tsv', comment = '#')
# featureCounts - feature
# tabfiles = c('../data/Pol2pausing/create_pausing_regions.test/featureCounts_dm3_ensembl78.pausing_tss.tsv', 
#              '../data/Pol2pausing/create_pausing_regions.test/featureCounts_dm3_ensembl78.pausing_tssds.tsv',
#              '../data/Pol2pausing/create_pausing_regions.test/featureCounts_dm3_ensembl78.pausing_txunit.tsv')
#
#
# for(tab in tabfiles){
#   cts.tab <- read_tsv(tab, comment = '#')
#   # featureCounts - gene
#   libsize0 <- colSums(genes.tab[colnames(cts.tab)[-c(1:6)]])
#   
#   cts.norm <- cbind(cts.tab[,1:6],tpm4(cts.tab[-c(1:6)], libsize = libsize0 , len = cts.tab$Length))
#   
#   write_tsv(x = cts.norm, gsub('.tsv$','.norm.tsv',tab))
# }
