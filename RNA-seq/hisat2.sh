#!/bin/bash
#SBATCH -J hisat2
#SBATCH -N 1 -c 80
#SBATCH -e /home/tmpdir/wangmingyang/RNAbam/log/job-%j_%a.err
#SBATCH -o /home/tmpdir/wangmingyang/RNAbam/log/job-%j_%a.log
#SBATCH -p fat
#SBATCH --mem=30Gb
#SBATCH -a 0-42
            


ulimit -n 65536
output=/home/tmpdir/wangmingyang/RNAbam/output
input=/home/tmpdir/wangmingyang/RNAbam/reads/clean
PBS_ARRAYID=${SLURM_ARRAY_TASK_ID} 
array_jobid=${SLURM_ARRAY_JOB_ID}  
jobid=${SLURM_JOB_ID}     
job_name=${SLURM_JOB_NAME} 
workdir=${SLURM_SUBMIT_DIR}
pid=${SLURM_TASK_PID}
nprocs=${SLURM_JOB_CPUS_PER_NODE}

uid="rna-wmy"
ls_date=`date +m%d%H%M%S`
temp=/tmpdisk/${uid}_${ls_date}_${PBS_ARRAYID}
mkdir ${temp}
cd ${temp}
pwd
bamfiles=(`ls ${input}/*_1.fq.gz`)
mybam=${bamfiles[$PBS_ARRAYID]}
sample=`basename ${mybam} _1.fq.gz`

hisat2 -p 80 -I 0 -X 500 --qc-filter -x $ref -1 ${input}/${sample}_1.fq.gz -2 ${input}/${sample}_2.fq.gz -S ${sample}.sam
samtools view -@ 80 -Sb ${sample}.sam > ${sample}.bam
rm -r ${sample}.sam
samtools sort -@ 80 -o ${sample}_sort.bam  ${sample}.bam
rm -r ${sample}.bam
samtools index -@ 80 ${sample}_sort.bam


mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

