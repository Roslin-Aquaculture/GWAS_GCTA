# GWAS using GCTA by Clemence Fraslin

Here we describe a step-by-step pipeline to perform a GWAS. The pipeline  uses the PLINK software for quality control and formatting of genotyping data, and the GCTA software to estimate genetic parameters and perform the GWAS , and. T the R package “QQMAN” is used to produce a Manhattan plot. The three software packages used in this pipeline are freely available and user-friendly. For more detailed information on the options available for each software, you can refer to the PLINK (http://zzz.bwh.harvard.edu/plink/ for version 1.07 or https://www.cog-genomics.org/plink/ for version 1.09) and GCTA websites (http://cnsgenomics.com/software/gcta/). The QQMAN package is available from CRAN (https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html) or GitHub (https://github.com/stephenturner/qqman).

## 0.	Optional: Merge datasets
If using genotyping performed on multiple SNP arrays (chip1 and chip2 here), you should merge PED files from all arrays with the “--merge” option:
```
plink --chr-set CHR --file mydata-chip1 --merge mydata-chip2.ped mydata.map --recode --out mydata
```
Before merging files, you should check if the map files are identical or use two different map files, and verify than alleles are coded equally (i.e., ATGC or 1234). If genotypes are coded differently, the options “--allele1234” or “--alleleACGT” have to be used before merging the datasets. 

## 1.	Format conversion
This first command transform .ped and .map data in binary files:
```
plink --ped mydata.ped --map mydata.map --chr-set CHR --make-bed --out newdata
```
You can add all the following options depending on your dataset “--allow-no-sex --nonfounders --no-fid --no-parents --no-sex --no-pheno --recode-allele refile.txt”.
The “--recode-allele” option is not mandatory but is useful to keep track of the reference alleles and make sure that alleles in future files are consistently coded. Indeed, PLINK recodes alleles according to the minor allele count (by default, the minor allele is A1 and the major allele is A2 in PLINK outputs), and as frequency of alleles might change between populations it is important to force a specific allele to ensure that the results are comparable. “refile.txt” is a two columns file were the first column is the SNP ID and the second column the allele to be used as reference allele (i.e. A1). The reference allele is typically the allele / nucleotide that appears in the reference genome sequence for the species of interest. 
Once your files are merged and/or in binary format you can start performing quality control. The command line to load binary files is --bfile.

## 2.	Removal of duplicated samples
This command is used to estimate pairwise IBD to detect duplicated samples:
```
plink --bfile newdata --genome --chr-set CHR --out tempdata
```
A file named tempdata.genome will be created with the IBD values of each pair of individuals in column 10. To remove duplicated samples, an extra step is necessary to select fish to be removed based on a threshold (IBD_thr = 0.90):
```
gawk -v IBD_thr=0.90 ' NR > 1 { if ($10 > IBD_thr) { print $1" " $2 " " $3 " " $4 }}' tempdata.genome > duplicated_ID.txt
gawk '{print $1 " " $2}' duplicated_ID.txt > ID_to_remove.txt
gawk '{print $3 " " $4}' duplicated_ID.txt >> ID_to_remove.txt
```

Those three lines of awk are used to create a file that contain the Family_ID and the fish_ID of duplicated samples to be removed from the dataset using the following plink line:
```
plink --bfile newdata --remove ID_to_remove.txt --make-bed --out noduplicdata
```

## 3.	Quality control
This command is used to perform all other quality controls on the dataset:
```
plink --bfile noduplicdata --chr-set CHR --geno 0.05 --mind 0.05 --hwe 0.000001 --maf 0.05 --nonfounders --make-bed --out finaldata
```

## 4.	GWAS
The following commands first create the GRM, perform the AI-REML analysis and run a GWAS using the MLMA approach with an example without covariate or fixed effect in the mode and an example with covariate and fixed effects:
```
gcta64 --bfile finaldata --autosome-num CHR --make-grm --out GRM_data --thread-num 10
gcta64 --reml --grm GRM_data --pheno TRAIT.pheno --qcovar FixedEffet.qcovar --covar Covariate.covar --out reml_TRAIT
gcta64 --mlma --bfile finaldata --grm GRM_data --autosome-num CHR --pheno TRAIT.pheno  --out gwas > gwas.log
gcta64 --mlma --bfile finaldata --grm GRM_data --pheno TRAIT.pheno --covar FixedEffet.covar --qcovar Covariate.qcovar --autosome-num CHR --out gwas-fixed-cov > gwas-fixed-cov.log
```
In the reml_TRAIT.hsq (output for --reml) file you will find estimates of h2, genetic and phenotypic variances. In the gwas.log file you can check that the REML analysis converged (typically AI-REML should converge in less than 20 iterations). 

If you want to perform a GWAS using the LOCO (leave-one-chromosome-out) approach, you can use the following command:
```
gcta64 --mlma-loco --bfile finaldata --grm GRM_data --pheno TRAIT.pheno --qcovar Covariate.qcovar --covar FixedEffect.covar --autosome-num CHR --out gwas-fixed-cov --thread-num 10 > gwas_loco-fixed-cov.log
```

## 5.	Manhattan plot
The following R commands use the QQMAN package to plot Manhattan plots from .mlma or .loco.mlma files:
```
install.package(“qqman”)
library(qqman)
NCHR=number_of_chromosomes   #Set "number_of_chromosomes" to the correct value for your species

#Load the GCTA output file
gwasResults = read.table(“gwas.mlma”, header = TRUE)

# Calculate the 5% significance level for the genome-wide ("Thrgen") and the chromosome-wide ("Thrchr") levels
Thrgen = -log10(0.05/nrow(gwasResults))
Thrchr = -log10(0.05/(nrow(gwasResults)/NCHR))

#Draw the Manhattan plot
manhattan(gwasResults, chr = "Chr", bp="bp", snp ="SNP", p="p", main = "Manhattan Plot", col = c("blue", "orange", “forestgreen”, “purple”), suggestiveline = Thrchr, genomewideline = Thrgen)
```
You can also use the qq function of the QQMAN package to plot the QQplot (quantile-quantile) of observed versus expected –log10(p-values):
```
qq(gwasResults$p, main = "Q-Q plot of GWAS p-values")
```

