input=/home/tmpdir/wangmingyang/ATAC
ref=/home/tmpdir/wangmingyang/ref
picard="singularity exec -B ${input}:${input} /home/tmpdir/wangmingyang/software/picard.sif"

$picard java -jar /usr/local/bin/picard/picard.jar AddOrReplaceReadGroups I=${input}/${sample}.bam O=${sample}_withRG.bam RGID=${sample} RGLB=Library1 RGPL=ILLUMINA RGPU=Unit1 RGSM=${sample} CREATE_INDEX=true

$picard java -jar /usr/local/bin/picard/picard.jar MarkDuplicates \
I=${sample}_withRG.bam O=${sample}_withRG_dedup.bam \
METRICS_FILE=${sample}_picard.dupMark.txt \
REMOVE_DUPLICATES=true CREATE_INDEX=true


sambamba="/home/tmpdir/wangmingyang/software/sambamba"
$sambamba markdup -r -t 10 ${sample}_withRG_dedup.bam ${sample}_withRG_dedup_pcr.bam
$sambamba index -t 10 ${sample}_withRG_dedup_pcr.bam


$samtools view -h ${sample}_withRG_dedup_pcr.bam | grep -v 'chrM' | samtools view -bS -o ${sample}_withRG_dedup_pcr_final.bam
$samtools index ${sample}_withRG_dedup_pcr_final.bam
