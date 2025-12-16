input=/home/tmpdir/wangmingyang/reads/single_cell
ref=/home/tmpdir/wangmingyang/ref
dnbc4tools="/home/wangmingyang/software/dnbc4tools/dnbc4tools"
data="/home/tmpdir/wangmingyang/single-cell/B001_cDNA"
data1="/home/tmpdir/wangmingyang/single-cell/B001_oligo"
$dnbc4tools rna run \
--name B001 \
--cDNAfastq1 ${data}/E250066871_L01_50_1.fq.gz,${data}/E250066871_L01_58_1.fq.gz,${data}/E250066871_L01_66_1.fq.gz \
--cDNAfastq2 ${data}/E250066871_L01_50_2.fq.gz,${data}/E250066871_L01_58_2.fq.gz,${data}/E250066871_L01_66_2.fq.gz \
--oligofastq1 ${data1}/B001_oligo-74_1.fq.gz \
--oligofastq2 ${data1}/B001_oligo-74_2.fq.gz \
--genomeDir ${input}/susScr11_pig \
--threads 30






