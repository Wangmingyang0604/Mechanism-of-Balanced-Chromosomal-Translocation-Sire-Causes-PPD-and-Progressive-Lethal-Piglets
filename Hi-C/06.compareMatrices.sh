#!/bin/bash
#SBATCH -J hic
#SBATCH -N 1 -c 10
#SBATCH -e /home/tmpdir/wangmingyang/HiC/log/job-%j_%a.err
#SBATCH -o /home/tmpdir/wangmingyang/HiC/log/job-%j_%a.log
#SBATCH -p fat
#SBATCH --mem=100Gb
#SBATCH -a 0
            


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
cooler="singularity exec  /home/tmpdir/wangmingyang/software/02_HiCtools.sif cooler"
hicConvertFormat="singularity exec -B /home/tmpdir/wangmingyang/HiC/D014output/results/hic_results/matrix/D014/raw:/mnt/wangmingyang /home/wangmingyang/software/hicexplorer_3_7_2.sif hicConvertFormat"
hicFindTADs="singularity exec -B /home/tmpdir/wangmingyang/HiC/:/home/tmpdir/wangmingyang/HiC/ /home/wangmingyang/software/hicexplorer_3_7_2.sif hicFindTADs"

#bamfiles=(`ls ${input}/*.h5`)
#mybam=${bamfiles[$PBS_ARRAYID]}
#sample=`basename ${mybam} .cool`


singularity exec -B /home/tmpdir/wangmingyang:/home/tmpdir/wangmingyang /home/wangmingyang/software/hicexplorer_3_7_2.sif hicCompareMatrices --matrices norm_chr9.h5 balance_chr9.h5 --outFileName diff_norm_balance_chr9.h5 --operation log2ratio 

mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

