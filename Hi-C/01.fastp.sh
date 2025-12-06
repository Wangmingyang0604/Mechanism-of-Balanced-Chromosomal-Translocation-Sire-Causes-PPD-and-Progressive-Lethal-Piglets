#!/bin/bash
#SBATCH -J hic
#SBATCH -N 1 -c 20
#SBATCH -e /home/tmpdir/Polytoe/HiC/log/job-%j_%a.err
#SBATCH -o /home/tmpdir/Polytoe/HiC/log/job-%j_%a.log
#SBATCH -p fat
#SBATCH --mem=30Gb
#SBATCH -a 0-2
            


ulimit -n 65536
output=/home/tmpdir/Polytoe/HiC/input
input=/home/tmpdir/reads  
PBS_ARRAYID=${SLURM_ARRAY_TASK_ID} 
array_jobid=${SLURM_ARRAY_JOB_ID}  
jobid=${SLURM_JOB_ID}     
job_name=${SLURM_JOB_NAME} 
workdir=${SLURM_SUBMIT_DIR}
pid=${SLURM_TASK_PID}
nprocs=${SLURM_JOB_CPUS_PER_NODE}

uid="hic-wmy"
ls_date=`date +m%d%H%M%S`
temp=/tmpdisk/${uid}_${ls_date}_${PBS_ARRAYID}
mkdir ${temp}
cd ${temp}
pwd
fastp="singularity exec /home/wangmingyang/software/fastp0.23.2.simg fastp"
bamfiles=(`ls ${input}/*_R1.fastq.gz`)
mybam=${bamfiles[$PBS_ARRAYID]}
sample=`basename ${mybam} _R1.fastq.gz`

$fastp -i  ${input}/${sample}_R1.fastq.gz -I ${input}/${sample}_R2.fastq.gz -o clean_${sample}_R1.fastq.gz -O clean_${sample}_R2.fastq.gz -h ${sample}_report.html -j ${sample}_report.json --thread 20


mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

