input=/home/wangmingyang/bamSus11.1/ONT
ref=/home/tmpdir/wangmingyang/ref/susScr11.fa

singularity exec -B ${input}:${input} \
/home/tmpdir/wangmingyang/software/figeno_FREEC_delly.sig delly lr -y ont -o ${sample}.bcf -g ${ref} ${input}/${sample}.bam
bcftools view ${sample}.bcf > ${sample}.vcf


singularity exec -B ${input}:${input} \
/home/tmpdir/wangmingyang/software/cuteSV.sif cuteSV \
-l 50 -s 10 -L 1000000 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
${input}/${sample}.bam ${ref} ${sample}.cuteSV.vcf ./

singularity exec -B ${input}:${input} \
/home/tmpdir/wangmingyang/software/sniffles.sif sniffles --input ${input}/${sample}.bam --vcf ${sample}_sniffles.vcf


/home/tmpdir/wangmingyang/software/SURVIVOR-1.0.6/Debug/SURVIVOR merge \
${input}/sample.txt 1000 2 1 1 0 100 ${sample}_merged.vcf

