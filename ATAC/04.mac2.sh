input=/home/tmpdir/wangmingyang/ATAC
ref=/home/tmpdir/wangmingyang/ref

$macs2 callpeak -t ${input}/${sample}.bam -n ${sample} --shift -100 --extsize 200 --nomodel -B --SPMR -g 2455392186 --outdir ${output} 2> ${sample}.macs2.log

