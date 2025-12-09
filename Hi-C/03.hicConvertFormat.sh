#!/bin/bash
#SBATCH -J hic
#SBATCH -N 1 -c 20
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
hicConvertFormat="singularity exec -B /home/tmpdir/wangmingyang/HiC/D014output/results/hic_results/matrix/D014/raw:/mnt/wangmingyang /home/wangmingyang/software/hicexplorer_3_7_2.sif hicConvertFormat"
#bamfiles=(`ls ${input}/*_R1.fastq.gz`)
#mybam=${bamfiles[$PBS_ARRAYID]}
#sample=`basename ${mybam} _R1.fastq.gz`


$hicConvertFormat  --matrices /mnt/wangmingyang/100000/D014_100000.matrix --bedFileHicpro /mnt/wangmingyang/100000/D014_1000000_abs.bed --inputFormat hicpro --outputFormat cool --outFileName D014_100000.cool --resolutions 100000

$hicConvertFormat  --matrices /mnt/wangmingyang/20000/D014_20000.matrix --bedFileHicpro /mnt/wangmingyang/20000/D014_20000_abs.bed --inputFormat hicpro --outputFormat cool --outFileName D014_20000.cool --resolutions 20000

$hicConvertFormat  --matrices /mnt/wangmingyang/40000/D014_40000.matrix --bedFileHicpro /mnt/wangmingyang/40000/D014_40000_abs.bed --inputFormat hicpro --outputFormat cool --outFileName D014_40000.cool --resolutions 40000


mkdir -p ${output}
cp -r ${temp}/* ${output}
echo "Remove tmp files at ${temp}"
rm -rf ${temp}

