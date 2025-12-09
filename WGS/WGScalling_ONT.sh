input=/home/tmpdir/wangmingyang/reads/ONT
ref=/home/tmpdir/wangmingyang/ref/susScr11.fa

RG="@RG\\tID:${sample}\\tPL:ONT\\tSM:${sample}\\tLB:${sample}\\tPU:1"
echo $RG
/home/jxlabgdp/01.Biosoft/bin/minimap2 -x map-hifi -R $RG -a -t 100 -y --secondary=no ${ref} ${input}/${sample}.fq.gz \
|samtools sort -@ 100  -o ${sample}_sort.bam

samtools index ${sample}_sort.bam

hifiasm="singularity exec -B ${input}:${input} /home/tmpdir/wangmingyang/software/hifiasm_0.25.0.sif hifiasm"

$hifiasm -t 80 --ont -o D018.asm ${input}/${sample}.fastq.gz




