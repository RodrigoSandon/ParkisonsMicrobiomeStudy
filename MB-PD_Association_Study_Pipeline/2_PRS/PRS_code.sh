### Polygenic risk score analyses of 119 bacteria genera versus PD risk

## Make a list of summary stats file 
ls *csv > /Volumes/T7Touch/NIHSummer2021/MB-PD_Association_Study_Pipeline/1_Pre-PRS/Part-2/genera.txt ## These files have been previously pruned by p-value < 0.05

## Format these files
cat genuses.txt | while read LINE 
do
  echo $LINE
  sed 's/\"//g' $LINE | sed 's/,/ /g' > temp.txt
  awk '{print $0"\t"$2":"$3}' temp.txt | sed 's/chr\:bp/ID/' > $LINE.temp_formatted.txt
  rm temp.txt
done

## Identify independent risk SNPs using our in-house LD reference data for European populations (/data/LNG/pdMeta5v2/cojoRef/binaryForCojo) - default clumping threshold (r2 = 0.1, kb= 250kb) + 10000 permutations

Rscript /Volumes/T7Touch/NIHSummer2021/MB-PD_Association_Study_Pipeline/2_PRS/PRSice.R --cov-file /Volumes/T7Touch/NIHSummer2021/MB-PD_Association_Study_Pipeline/2_PRS/sample_covariate.txt -t /data/LNG/saraB/concat_HARDCALLS_PD_september_2018_no_cousins --beta --snp ID --A1 eff.allele --A2 ref.allele --stat beta --se SE --pvalue pDerived --ld /data/LNG/pdMeta5v2/cojoRef/binaryForCojo  --print-snp --score std --perm 10000 --prsice /data/LNG/pdMeta5v2/leaveOneOutPrsice/PRSice_linux/PRSice_linux -n 8 --binary-target T --quantile 4 --prevalence 0.005 --fastscore -b /Volumes/T7Touch/NIHSummer2021/MB-PD_Association_Study_Pipeline/1_Pre-PRS/Part-1/sumstats-p/addedP_bac_1.csv.temp_formatted.txt --out bac_1_prs

## We are only using two of the 119 bacterial genera (our PRS hits) in this example to avoid redundancy.

## Remove NeuroX individuals & extract nominated variants
cat genuses_formatted_list.txt | while read LINE 
do
plink --bfile /data/LNG/saraB/concat_HARDCALLS_PD_september_2018_no_cousins --remove-fam NeuroX.txt --extract $LINE.snps --make-bed --out pruned_$LINE
done

## Make score files
cat genuses_formatted_list.txt | while read LINE 
do
awk '{print $14, $6, $7}' addedPgenus.$LINE.summary.txt.csv.temp_formatted.txt | sed '1d' > $LINE.toscore.txt
done

## Make sure score files have 3 expected fields rather than 2
cat genuses_formatted_list.txt | while read LINE 
do
grep ":" $LINE.toscore.txt > true_$LINE.toscore.txt
done

##awk: fatal: cannot open file `addedPgenus.Clostridiuminnocuumgroup.id.14397.summary.txt.csv.temp_formatted.txt' for reading (No such file or directory)
##awk: fatal: cannot open file `addedPgenus..summary.txt.csv.temp_formatted.txt' for reading (No such file or directory)

## Calculate scores
cat genuses_formatted_list.txt | while read LINE 
do
plink --bfile pruned_$LINE --score $LINE.toscore.txt --make-bed --out pruned_$LINE
done

## Run PRS (logistic regression) on R
R
library("data.table")
listOfProfiles <- read.table("genuses_formatted_list.txt", header = T)
names(listOfProfiles) <- c("id")
covs1 <- fread("/data/LNG/saraB/WGS/noage_toPRSice_phenosAndCovs_renamed.tab", header = T)
covs2 <- fread("/data/LNG/saraB/concat_HARDCALLS_PD_september_2018_no_cousins.fam", header = F)
colnames(covs2) <- c("FID", "IID", "MAT", "PAT", "SEX", "PHENO")
covsfinal <- merge (covs1, covs2, by ="FID")
covsfinal$CASE <- covsfinal$PHENO.x - 1
outPut <- matrix(ncol = 4, nrow = length(listOfProfiles$id), NA)
colnames(outPut) <- c("genus","b","se","p")
for(i in 1:length(listOfProfiles$id))
{
	profileName <- 	as.character(listOfProfiles$id[i])
	profile <- fread(file = paste(profileName, ".profile", sep = ""), header = T)
	profile$index <- paste(profile$FID, profile$IID, sep = "")
	data <- merge(covsfinal, profile, by = "index")
	meanControls <- mean(data$SCORE[data$CASE == 0])
	sdControls <- sd(data$SCORE[data$CASE == 0])
	data$zSCORE <- (data$SCORE - meanControls)/sdControls	
	grsTest <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + DUTCH + FINLAND + GERMANY + HBS + MCGILL + MF + NIA + OSLO + PDBP + PPMI + SHULMAN + SPAIN3 + TUBI, family="binomial", data = data)	
	beta <- summary(grsTest)$coefficients["zSCORE","Estimate"]
	se <- summary(grsTest)$coefficients["zSCORE","Std. Error"]
	p <- summary(grsTest)$coefficients["zSCORE","Pr(>|z|)"]
	outPut[i,1] <- profileName
	outPut[i,2] <- beta
	outPut[i,3] <- se
	outPut[i,4] <- p
}
write.table(outPut, "Genus_PRS.tab", quote = F, sep = "\t", row.names = F)

## TEST ANALYSIS for Actinomyces

# library("data.table")
# temp_data <- read.table("pruned_Actinomyces.id.423.profile", header = T) 
# covs1 <- fread("/data/LNG/saraB/WGS/noage_toPRSice_phenosAndCovs_renamed.tab", header = T)
# covs2 <- fread("/data/LNG/saraB/concat_HARDCALLS_PD_september_2018_no_cousins.fam", header = F)
# colnames(covs2) <- c("FID", "IID", "MAT", "PAT", "SEX", "PHENO")
# covsfinal <- merge (covs1, covs2, by ="FID")
# data <- merge(temp_data, covsfinal, by = "FID")
# data$CASE <- data$PHENO.x - 1
# meanControls <- mean(data$SCORE[data$CASE == 0])
# sdControls <- sd(data$SCORE[data$CASE == 0])
# data$zSCORE <- (data$SCORE - meanControls)/sdControls
# grsTests <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = data)
# summary(grsTests)
# 
# 
# Call:
# glm(formula = CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + 
#     PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial", data = data)
# 
# Deviance Residuals: 
#    Min      1Q  Median      3Q     Max  
# -1.919  -1.003  -0.810   1.278   1.796  
# 
# Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.43040    0.03974  10.831  < 2e-16 ***
# zSCORE        0.03021    0.01424   2.121   0.0339 *  
# SEX          -0.58575    0.02603 -22.505  < 2e-16 ***
# PC1         -35.38600    2.73911 -12.919  < 2e-16 ***
# PC2          50.17593    2.86877  17.490  < 2e-16 ***
# PC3          10.63757    2.68575   3.961 7.47e-05 ***
# PC4           0.63048    2.65991   0.237   0.8126    
# PC5          16.19265    2.70187   5.993 2.06e-09 ***
# PC6         -24.20283    2.74525  -8.816  < 2e-16 ***
# PC7           1.61607    2.62627   0.615   0.5383    
# PC8          12.66905    2.70987   4.675 2.94e-06 ***
# PC9          -5.96044    2.64009  -2.258   0.0240 *  
# PC10        -16.18338    2.65955  -6.085 1.16e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 35275  on 26385  degrees of freedom
# Residual deviance: 34124  on 26373  degrees of freedom
# AIC: 34150
# 
# Number of Fisher Scoring iterations: 4