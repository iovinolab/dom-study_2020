set -e
module load slurm umi_tools/1.0.0 cutadapt/2.5 bowtie2/2.3.3.1 sambamba/0.7.0 samtools/1.10 deeptools/3.3.1 MultiQC/1.8 seqtk/1.2 FastQC/0.11.5 bedops/2.4.27 bedtools2/2.27.0 subread/1.5.3 R/3.5.0

read directory config params <<< "$@"

mkdir -p $directory/cluster_log
snakemake -d $directory --configfile $config $params # --cluster "SlurmEasy -n {rule} -t {threads} -l $directory/cluster_log" --latency-wait 500 --keep-going
