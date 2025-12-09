input=/home/tmpdir/wangmingyang/illumina_cnv
ref=/home/tmpdir/wangmingyang/ref/susScr11.fa

/home/tmpdir/wangmingyang/FREEC-11.6b/src/freec -conf  ${input}/freec_WGS_control.config

singularity exec /home/tmpdir/wangmingyang/software/cnvkit_latest.sif cnvkit.py access ${ref} -s 10000 -o access-10kb.susScr11.bed

singularity exec -B ${input}:${input} \
/home/tmpdir/wangmingyang/software/cnvkit_latest.sif cnvkit.py batch \
/mnt/zhiyan/Polydactyly1.bam --normal ${input}/norm1.bam  -m wgs --fasta ${ref} \
--annotate ${input}/refFlat.txt \
--access ${input}/access-10kb.susScr11.bed \
--output-reference reference.cnn --output-dir ./


