geno_raw <- read.table("gene.txt", header=F)
colnames(geno_raw) <- c("SNP1", "SNP2", "SNP3") 
table(geno_raw$SNP1) 
table(geno_raw$SNP2)
table(geno_raw$SNP3)
pheno <- scan("phe.txt")  
if(nrow(geno_raw) != length(pheno)) stop("no match")
model <- lm(pheno ~ SNP1 + SNP2 + SNP3, data=geno_raw)
summary(model)
residuals <- model$residuals
mean_phenotype <- mean(pheno)

new_phenotype <- residuals + mean_phenotype

print(new_phenotype)

write.table(new_phenotype, file = "new_phenotype.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
