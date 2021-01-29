suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-t", "--table"), default=NA, type='character',
              help="featureCounts table (tsv)"),
  make_option(c("-s", "--samplesheet"), default=NA, type='character',
              help="Sample sheet defining \'cond\' and \'ctrl\' (tsv)"),
  make_option(c("-o", "--outprefix"), default=NA, type='character',
              help="output prefix"),

  make_option("--gtf", default=NA, type='character',
              help="Optional annotation file to map gene_id to symbol (default: %default)"),
  make_option("--gene_id", default='gene_id', type='character',
              help="Identifier for \'gene_id'\ in gtf (default: %default)"),
  make_option("--gene_name", default='gene_name', type='character',
              help="Identifier for \'gene_name'\ in gtf (default: %default)"),
  make_option("--batch", action = 'store_true', default=FALSE, type='logical',
              help="Use \'batch\' as covariate"),
  make_option("--scaling_geneset", default = NULL, type='character',
              help="Table with gene_id. Column name = 'gene_id', tsv (default: %default)"),
  

  make_option(c("--min_reads"), default = 10, type='integer',
              help="Discard genes with less than n reads on average (n = %default)"),

  make_option(c("--lfcShrinkage_method"), default = 'normal', type='character',
              help="DESeq2::lfcShrinkage method, must of \'normal\' or \'apeglm\' (default: %default)"),
  
  make_option(c("--skip_plots"), action = 'store_true', default = FALSE, type='logical',
              help="skip QC plots"),

  make_option(c("--lib.location"), default = NULL, type = 'character',
              help="Specify library location")
)
opt = parse_args(OptionParser(option_list=option_list))


lib.loc = opt$lib.location

suppressPackageStartupMessages(library(DESeq2, lib.loc = lib.loc))
suppressPackageStartupMessages(library(readr, lib.loc = lib.loc))
suppressPackageStartupMessages(library(janitor, lib.loc = lib.loc))
suppressPackageStartupMessages(library(reshape2, lib.loc = lib.loc))
suppressPackageStartupMessages(library(ggplot2, lib.loc = lib.loc))
suppressPackageStartupMessages(library(dplyr, lib.loc = lib.loc))

# adopted from
# https://www.biostars.org/p/335187/
tpm3_stable <- function(counts,len) {
  x <- log(counts) - log(len)
  return(exp(t(t(x) - log(colSums(exp(x))) + log(1e6))))
}

samplesheet <- read_tsv(opt$samplesheet, comment = '#')
colDat0 <- samplesheet
sample0 = make_clean_names(colDat0$sample)
samplename.map = setNames(colDat0$sample, nm = sample0)
colDat0$sample <- sample0
colDat0$condition <- relevel(as.factor(colDat0$condition), 'ctrl')

cts.tab <- read_tsv(opt$table, comment = '#')
colnames(cts.tab) <- make_clean_names(gsub('(.*)\\.bam','\\1',basename(colnames(cts.tab))))

cts.mat <- cts.tab[,-c(1:6)]
#
# Writing the design, according to DESeq2 vignettes
#
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Note: In order to benefit from the default settings of the package, you should
# put _the variable of interest at the end of the formula_ and make sure the
# control level is the first level.
if(opt$batch){
  dds <- DESeq2::DESeqDataSetFromMatrix(cts.mat[colDat0$sample], colData = colDat0, design = ~ batch + condition)
} else {
  dds <- DESeq2::DESeqDataSetFromMatrix(cts.mat[colDat0$sample], colData = colDat0, design = ~ condition)
}

rownames(dds) <- cts.tab$geneid
mcols(dds)$basepairs <- cts.tab$length

# filter by read counts
keep <- rowMeans(assay(dds)) > opt$min_reads
dds <- dds[keep]

if(!is.null(opt$scaling_geneset)){
  cat('>> Computing sizeFactors from gene set\n\t',opt$scaling_geneset,'\n')
  scaling_tab = read_tsv(opt$scaling_geneset)
  gene_id.scaling = intersect(scaling_tab$gene_id, rownames(dds))
  
  dds0 = DESeq2::estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_id.scaling)
  
  sizeFactors(dds) = sizeFactors(dds0)
}

dds <- DESeq(dds)
coef0 = grep('^condition_', resultsNames(dds), value = TRUE)
cat("Report for coefficient:\n", coef0,'\n')
res <- data.frame(gene_id = rownames(dds),
                  results(dds, name = coef0)) %>% 
  mutate(significant = ifelse(padj < 0.05, ifelse(log2FoldChange > 0, 'up','down') , 'none'))

res.lfc_shrinkage = lfcShrink(dds, coef = coef0, 
                              type = opt$lfcShrinkage_method) %>% 
  as.data.frame()
res.lfc_shrinkage = data.frame(gene_id = rownames(res.lfc_shrinkage), res.lfc_shrinkage, row.names = NULL) %>% 
  mutate(significant = ifelse(padj < 0.05, ifelse(log2FoldChange > 0, 'up','down') , 'none'))

cat('> Computing rlog and rpkm\n')
rlog0 <- rlog(dds, blind = FALSE)
rpkm0 <- fpkm(dds)

if(!is.na(opt$gtf)){
  cat('> Adding symbols\n')
  suppressPackageStartupMessages(library(rtracklayer, lib.loc = lib.loc))
  suppressPackageStartupMessages(library(tibble, lib.loc = lib.loc))
  anno.gtf <- import.gff(opt$gtf)
  gid2symbol <- deframe(mcols(anno.gtf)[c(opt$gene_id, opt$gene_name)])
  
  cat('>> to DESeq2 results\n')
  res <- res %>% mutate(symbol = gid2symbol[as.character(gene_id)])
  res.lfc_shrinkage = res.lfc_shrinkage %>% mutate(symbol = gid2symbol[as.character(gene_id)])
  
  cat('>> to RPKM table\n')
  rpkm0 <- rpkm0 %>% as.data.frame %>%
    mutate(gene_id = as.character(rownames(rpkm0)),
           symbol = gid2symbol[as.character(rownames(rpkm0))])
  b = colnames(rpkm0) %in% names(samplename.map)
  colnames(rpkm0)[b] = samplename.map[colnames(rpkm0)[b]]
}


# cat(">>> Colnames\n")
# print(colnames(res))
# print(colnames(rpkm0))
c.cols = c('gene_id','symbol')

write_rds(x = list(dds = dds, results = res, rlog = rlog0, rpkm = rpkm0), path = paste0(opt$outprefix, '.DESeq2.rds'))

write_tsv(x = res[c(c.cols[1], setdiff(colnames(res),c.cols), c.cols[2])], paste0(opt$outprefix,'.DESeq2.tsv'))

write_tsv(x = rpkm0[c(c.cols, setdiff(colnames(rpkm0),c.cols))], paste0(opt$outprefix,'.rpkm_table.tsv'))

d.sizeFactors = data.frame(sample = samplename.map[names(sizeFactors(dds))], sizeFactor = sizeFactors(dds), row.names = NULL)
write_tsv(x = d.sizeFactors, paste0(opt$outprefix,'.sizeFactors.tsv'))

# cat('---\n')
if(!opt$skip_plots){
  cat('> Computing QC plots\n')
  cat('>> Computing MA plot\n')
  p.ma <- ggplot(res %>% arrange(padj < 0.05), aes(log2(baseMean), log2FoldChange)) + geom_point(aes(color = significant), show.legend = FALSE) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c('up' = 'red','down' = 'red','none' = 'grey'))
  ggsave(paste0(opt$outprefix,'.ma_plot.png'), p.ma)

  cat('>> Computing Volcano plot\n')
  p.volcano <- ggplot(res %>% arrange(padj < 0.05), aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(color = significant), show.legend = FALSE) +
    geom_vline(xintercept = 0) +
    scale_color_manual(values = c('up' = 'red','down' = 'red','none' = 'grey'))
  ggsave(paste0(opt$outprefix,'.volcano_plot.png'), p.volcano)

  cat('>> Computing PCA plot\n')
  p.pca <- plotPCA(object = rlog0) + ggtitle('GRO-seq PCA - rlog normalized counts')
  p.pca_data <- plotPCA(object = rlog0, returnData = TRUE)
  p.pca_data$name = samplename.map[p.pca_data$name]
  
  ggsave(paste0(opt$outprefix,'.pca_plot.png'), p.pca)
  write_tsv(p.pca_data, paste0(opt$outprefix,'.pca_plot.data.tsv'))

  cat('>> Computing sample distance plot\n')
  sampleDists <- dist(t(assay(rlog0)))
  library("RColorBrewer")
  library("pheatmap")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(rlog0)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           display_numbers = TRUE,
           filename = paste0(opt$outprefix, '.sample_distance.png'))
}

graphics.off()

cat('Output written to', opt$outprefix ,'\n')