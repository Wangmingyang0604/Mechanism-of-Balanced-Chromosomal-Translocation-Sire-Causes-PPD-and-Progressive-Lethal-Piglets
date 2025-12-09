#!/bin/bash
#SBATCH -J stringtie
#SBATCH -N 1 -c 80
#SBATCH -e /home/tmpdir/wangmingyang/RNAbam/log/job-%j_%a.err
#SBATCH -o /home/tmpdir/wangmingyang/RNAbam/log/job-%j_%a.log
#SBATCH -p fat
#SBATCH --mem=30Gb
#SBATCH -a 0-42
            


ulimit -n 65536
output=/home/tmpdir/wangmingyang/RNAbam/bam/output
input=/home/tmpdir/wangmingyang/RNAbam/bam
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
bamfiles=(`ls ${input}/*.bam`)
mybam=${bamfiles[$PBS_ARRAYID]}
sample=`basename ${mybam} .bam`

bedtools coverage -a ${input}/STK17A.bed -b ${input}/${sample}.bam -counts > ${sample}_STK17A_coverage.bed

bedtools coverage -a ${input}/TYR.bed -b ${input}/${sample}.bam -counts > ${sample}_TYR_coverage.bed

mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

