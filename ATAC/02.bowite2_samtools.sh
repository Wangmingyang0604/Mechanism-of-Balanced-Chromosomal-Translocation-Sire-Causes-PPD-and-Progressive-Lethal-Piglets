input=/home/tmpdir/wangmingyang/ATAC
ref=/home/tmpdir/wangmingyang/ref
samtools="/hpai/aios3.0/software/anaconda3/bin/samtools"
bowtie2="singularity exec -B ${input}:${input} /home/wangmingyang/software/hicpro.sif bowtie2"

$bowtie2 --very-sensitive -X 2000 -p 40 -x ${ref}/susScr11  -1 ${input}/${sample}_1.fastq.gz   -2 ${input}/${sample}_2.fastq.gz -S ${sample}.sam 2> ${sample}_bowtie2.txt



$samtools view -@ 40 -bS ${input}/${sample}.sam | samtools sort -@ 40 -o ${sample}_sort.bam
$samtools index -@ 40 ${sample}_sort.bam
$samtools view -@ 40 -b -F 524 -q 30 ${sample}_sort.bam -o filtered_${sample}_sort.bam
$samtools view -@ 40 -c filtered_${sample}_sort.bam

