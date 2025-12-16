input=/home/tmpdir/wangmingyang/reads/single_cell
ref=/home/tmpdir/wangmingyang/ref
cellbender="singularity exec -B /home/tmpdir/wangmingyang:/home/tmpdir/wangmingyang /home/tmpdir/wangmingyang/software/cellbender
data="/home/tmpdir/wangmingyang/single-cell/prepare/matrix_result/B001/outs/raw_matrix/huada"
$cellbender cellbender remove-background \
--input ${data}/ \
--output B001_cell_bender.h5 \
--cpu-threads 28





