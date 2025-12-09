input=/home/tmpdir/wangmingyang/reads/DNA
ref=/home/tmpdir/wangmingyang/ref/susScr11.fa

fastp="singularity exec /home/wangmingyang/software/fastp0.23.2.simg fastp"
bamfiles=(`ls ${input}/*_R1.fastq.gz`)
mybam=${bamfiles[$PBS_ARRAYID]}
sample=`basename ${mybam} _R1.fastq.gz`

$fastp -i  ${input}/${sample}_R1.fastq.gz -I ${input}/${sample}_R2.fastq.gz -o clean_${sample}_R1.fastq.gz -O clean_${sample}_R2.fastq.gz -h ${sample}_report.html -j ${sample}_report.json --thread 20

###mapping
bwa mem -t 50 ${ref} clean_${sample}_R1.fastq.gz clean_${sample}_R2.fastq.gz |samtools view -Sb > ${sample}.bam

#sort
samtools sort -@ 20 -T ${sample}.bam -m2G -O bam -o ${sample}.sorted.bam  ${sample}.bam

gatk AddOrReplaceReadGroups -I ${sample}.sorted.bam -O ${sample}.sorted.add.bam -ID $sample -LB library1 -PL illumina -SM $sample -PU unit1 &&
gatk MarkDuplicates -I ${sample}.sorted.add.bam -O ${sample}.sorted.add.markdup.bam -M ${sample}.sorted.add.markdup_metrics.txt &&
samtools index ${sample}.sorted.add.markdup.bam

chroms=($(awk '{print $1}' susScr11.fa.fai |sed 's/>//' |grep "chr"  | tr '\n' ' ' ))

for i in ${chroms[@]}
do
	$gatk HaplotypeCaller --java-options -Xmx10G  --native-pair-hmm-threads 10 -R $ref -I ${sample}.sorted.add.markdup.bam  -O ${sample}.sorted.add.${i}.markdup.g.vcf  -ERC GVCF -L ${i} &
done && wait

gvcfs=/home/tmpdir/wangmingyang/reads/DNA/input.list
outdir=/home/tmpdir/wangmingyang/reads/DNA/gatk_GenomicsDBImport
tmpdir=/home/tmpdir/wangmingyang/reads/DNA/tmp
outname="merge"

for i in ${chroms[@]}; do
    gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport \
        -R $ref \
        -V $gvcfs \
        -L $i \
        --genomicsdb-workspace-path $outdir/database_$i \
        --tmp-dir $tmpdir \
        --reader-threads 3 \
        --batch-size 50 &
done

for i in ${chroms[@]}; do
    gatk GenotypeGVCFs \
        -R $ref \
        -V gendb://$outdir/database_$i \
        -O $outdir/${outname}.$i.vcf &
done


bgzip merge.${chr}.vcf
tabix -p vcf merge.${chr}.vcf.gz
gatk SelectVariants -V merge.${chr}.vcf.gz -O merge.${chr}.snp.vcf --select-type-to-include SNP
gatk VariantFiltration -O merge.${chr}.vcf.temp -V merge.${chr}.snp.vcf \
--filter-expression 'QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
--filter-name lowQualFilter --cluster-window-size 10 --cluster-size 3 \
--missing-values-evaluate-as-failing







