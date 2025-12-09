plink --bfile poly --geno 0.05 --make-bed --out poly2
plink --bfile poly2 --maf 0.01 --make-bed --out poly3
plink --bfile poly3 --hwe 1e-6 --make-bed --out poly4

gemma -bfile poly4 -p phe.txt -gk 2
gemma -bfile Poly4 -k ${input}/output/result.sXX.txt -c cov.txt -lmm 2 -p phe.txt -o gwas_result


####condition-GWAS
plink --bfile Poly4 --snps chr18:51384660:IG --recode12 --out test
awk '{print 1,$7,$8}' test.ped >1.txt
awk '{print ($1==1 && $2==2 && $3==2) ? "1 0 0" : ($1==1 && $2==1 && $3==2) ? "0 1 0" : ($1==1 && $2==1 && $3==1) ? "0 0 1" : $0}' 1.txt > gene.txt
gemma --bfile Poly4 -k ${input}/output/result.sXX.txt -c cov2.txt -lmm 2 -p new_phenotype.txt -o condition_result
