#!/bin/bash
#SBATCH -J hic
#SBATCH -N 1 -c 20
#SBATCH -e /home/tmpdir/wangmingyang/HiC/log/job-%j_%a.err
#SBATCH -o /home/tmpdir/wangmingyang/HiC/log/job-%j_%a.log
#SBATCH -p fat
#SBATCH --mem=200Gb
#SBATCH -a 0-2
            


ulimit -n 65536
output=/home/tmpdir/wangmingyang/HiC/input
input=/home/tmpdir/wangmingyang/reads/clean
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
hicpro="singularity exec /home/wangmingyang/software/hicpro.sif /HiC-Pro_3.1.0/bin/HiC-Pro"
bamfiles=(`ls ${input}/*_R1.fastq.gz`)
mybam=${bamfiles[$PBS_ARRAYID]}
sample=`basename ${mybam} _R1.fastq.gz`

echo "HiC-Pro starting..."
mkdir data
mkdir -p data/hicpro
mkdir Scripts
mkdir Ref
cp ${input}/config-hicpro.txt ./Scripts
cp -r /home/wangmingyang/contain-Reference/ref ./Ref
ln -s ${input}/${sample}_R1.fastq.gz  ./data/hicpro/${sample}_R1.fastq.gz
ln -s ${input}/${sample}_R2.fastq.gz  ./data/hicpro/${sample}_R2.fastq.gz

$hicpro  -i data/  -o results/ -c Scripts/config-hicpro.txt &&


mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

