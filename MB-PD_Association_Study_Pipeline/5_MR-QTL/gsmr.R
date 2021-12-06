library("gsmr")
cand_sum_stats <- read.csv("addedPgenus.CandidatusSoleaferrea.id.11350.summary.txt.csv", header=T, sep=",")
head(cand_sum_stats)
dim(cand_sum_stats)

# SNP: the genetic instrument
# a1: effect allele
# a2: the other allele
# a1_freq: frequency of a1
# bzx: the effect size of a1 on risk factor
# bzx_se: standard error of bzx
# bzx_pval: p value for bzx
# bzx_n: per-SNP sample size of GWAS for the risk factor
# bzy: the effect size of a1 on disease  - DONT HAVE
# bzy_se: standard error of bzy  - DONT HAVE
# bzy_pval: p value for bzy  - DONT HAVE
# bzy_n: per-SNP sample size of GWAS for the disease