input=/home/tmpdir/wangmingyang/reads/single_cell
ref=/home/tmpdir/wangmingyang/ref
dnbc4tools="/home/wangmingyang/software/dnbc4tools/dnbc4tools"
$dnbc4tools tools mkgtf --action stat --ingtf ${ref}/Sus_scrofa.Sscrofa11.1.113.gtf --output gtfstat.txt --type gene_biotype

$dnbc4tools tools mkgtf --action check --ingtf ${ref}/Sus_scrofa.Sscrofa11.1.113.gtf --output corrected.gtf

$dnbc4tools tools mkgtf --ingtf corrected.gtf --output genes.filter.gtf --type gene_biotype




