module load R/3.5.2 subread/2.0.0 
set -e

cdir=$(dirname $0)

read nproc outprefix annotation_gtf bams <<< "$@"

# config needs
# gtf
# annotation_gtf=$(grep 'annotation_gtf:' $config_yaml | cut -f 2 -d ':')

segmentation_yaml=${outprefix}.yaml
Rscript $cdir/Segmentation.R --gtf $annotation_gtf -o $segmentation_yaml &> ${outprefix}.segmentation.log


declare -A regions
regions["tss"]="$(grep 'tss:' $segmentation_yaml | cut -f 2 -d ':')"
regions["tssds"]="$(grep 'tssds:' $segmentation_yaml | cut -f 2 -d ':')"
regions["txunit"]="$(grep 'txunit:' $segmentation_yaml | cut -f 2 -d ':')"
regions["gene_body"]="$(grep 'gene_body:' $segmentation_yaml | cut -f 2 -d ':')"

params="-T $nproc -Q 3 -M -s 1 -O -t transcript -g transcript_id"
for r0 in "tss" "tssds" "txunit";
do
  featureCounts $params -a ${regions[$r0]} -o ${outprefix}.${r0}.counts.tsv $bams &> ${outprefix}.${r0}.log
done

params_genes="-T $nproc -t gene -g gene_id -Q 3 -M -s 1 -O"
featureCounts $params_genes -a ${regions['gene_body']} -o ${outprefix}.genebody.counts.tsv $bams &> ${outprefix}.genebody.log


for r0 in "tss" "tssds" "txunit";
do
  Rscript $cdir/GROseq_normalization.R -c ${outprefix}.${r0}.counts.tsv -t ${outprefix}.genebody.counts.tsv -o ${outprefix}.${r0}.normalized.tsv &> ${outprefix}.${r0}.normalize.log
done

Rscript $cdir/GROseq_pausing.R -g $annotation_gtf -p ${outprefix}.tss.normalized.tsv -d ${outprefix}.txunit.normalized.tsv --outprefix ${outprefix}.bg_txunit.Pausing
Rscript $cdir/GROseq_pausing.R -g $annotation_gtf -p ${outprefix}.tss.normalized.tsv -d ${outprefix}.tssds.normalized.tsv --outprefix ${outprefix}.bg_tssds.Pausing

touch ${outprefix}.DONE
