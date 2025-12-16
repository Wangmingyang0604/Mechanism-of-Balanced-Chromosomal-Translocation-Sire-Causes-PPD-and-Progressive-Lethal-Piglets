#!/bin/bash
#SBATCH -J atac
#SBATCH -N 1 -c 20
#SBATCH -e /home/tmpdir/wangmingyang/ATAC/log/job-%j_%a.err
#SBATCH -o /home/tmpdir/wangmingyang/ATAC/log/job-%j_%a.log
#SBATCH -p fat
#SBATCH --mem=30Gb
#SBATCH -a 0-1
            


ulimit -n 65536
output=/home/tmpdir/wangmingyang/ATAC/output
input=/home/tmpdir/wangmingyang/ATAC
PBS_ARRAYID=${SLURM_ARRAY_TASK_ID} 
array_jobid=${SLURM_ARRAY_JOB_ID}  
jobid=${SLURM_JOB_ID}     
job_name=${SLURM_JOB_NAME} 
workdir=${SLURM_SUBMIT_DIR}
pid=${SLURM_TASK_PID}
nprocs=${SLURM_JOB_CPUS_PER_NODE}

uid="atac-wmy"
ls_date=`date +m%d%H%M%S`
temp=/tmpdisk/${uid}_${ls_date}_${PBS_ARRAYID}
mkdir ${temp}
cd ${temp}
pwd
fastp="singularity exec /home/wangmingyang/software/fastp0.23.2.simg fastp"
bamfiles=(`ls ${input}/*_1.fastq.gz`)
mybam=${bamfiles[$PBS_ARRAYID]}
sample=`basename ${mybam} _1.fastq.gz`

$fastp -i  ${input}/${sample}_1.fastq.gz -I ${input}/${sample}_2.fastq.gz -o clean_${sample}_1.fastq.gz -O clean_${sample}_2.fastq.gz -h ${sample}_report.html -j ${sample}_report.json --thread 20


mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

