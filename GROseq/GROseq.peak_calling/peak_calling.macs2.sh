# Requires the fragments in the bam files to be split by strand

module load MACS2/2.1.2 sambamba/0.7.0 bedtools2/2.27.0
nproc=8

outdir="MACS2/" && mkdir -p $outdir

bams_stranded="bams_stranded"
sample_wt_fwd=$(ls $bams_stranded/*T180*.fwd.bam)
sample_wt_rev=$(ls $bams_stranded/*T180*.rev.bam)

sample_cond_fwd=$(ls $bams_stranded/*T216*.fwd.bam)
sample_cond_rev=$(ls $bams_stranded/*T216*.rev.bam)

params="--nomodel --extsize 100 -q 0.05 --call-summits"

# forward strand
macs2 callpeak -n GRO-seq.wildtype.fwd $params -t $sample_wt_fwd -g dm --outdir $outdir/peaks.wildtype/stranded &> $outdir/GRO-seq.wildtype.fwd.log
# reverse strand
macs2 callpeak -n GRO-seq.wildtype.rev $params -t $sample_wt_rev -g dm --outdir $outdir/peaks.wildtype/stranded &> $outdir/GRO-seq.wildtype.rev.log

tmp_dir="$outdir/peaks.wildtype_tmp" && mkdir -p $tmp_dir
awk -v OFS='\t' '{$6="+"; print($0)}' $outdir/peaks.wildtype/stranded/GRO-seq.wildtype.fwd_peaks.narrowPeak  > $tmp_dir/peaks.fwd.narrowPeak &&
awk -v OFS='\t' '{$6="-"; print($0)}' $outdir/peaks.wildtype/stranded/GRO-seq.wildtype.rev_peaks.narrowPeak > $tmp_dir/peaks.rev.narrowPeak &&
cat $tmp_dir/peaks.fwd.narrowPeak $tmp_dir/peaks.rev.narrowPeak | sortBed > "./$outdir/peaks.wildtype.stranded.narrowPeak"


# forward strand
macs2 callpeak -n GRO-seq.condition.fwd $params -t $sample_cond_fwd -g dm --outdir $outdir/peaks.condition/stranded &> $outdir/GRO-seq.condition.fwd.log
# reverse strand
macs2 callpeak -n GRO-seq.condition.rev $params -t $sample_cond_rev -g dm --outdir $outdir/peaks.condition/stranded &> $outdir/GRO-seq.condition.rev.log

tmp_dir="$outdir/peaks.condition_tmp" && mkdir -p $tmp_dir
awk -v OFS='\t' '{$6="+"; print($0)}' $outdir/peaks.condition/stranded/GRO-seq.condition.fwd_peaks.narrowPeak  > $tmp_dir/peaks.fwd.narrowPeak &&
awk -v OFS='\t' '{$6="-"; print($0)}' $outdir/peaks.condition/stranded/GRO-seq.condition.rev_peaks.narrowPeak > $tmp_dir/peaks.rev.narrowPeak &&
cat $tmp_dir/peaks.fwd.narrowPeak $tmp_dir/peaks.rev.narrowPeak | sortBed > "./$outdir/peaks.condition.stranded.narrowPeak"
