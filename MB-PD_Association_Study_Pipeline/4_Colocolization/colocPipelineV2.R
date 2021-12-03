library("data.table")
install.packages("coloc") 
source("https://bioconductor.org/biocLite.R")
biocLite("snpStats")
install.packages("robustbase")
library("robustbase")
library("coloc")
library("tidyverse")

#The file that contains SNPs with significance of 5E-N for both disease and phenotype
masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/p_smr_multi_mrQTL_belowFDR"
files <- list.files(path=masterDir, pattern="*.tab", full.names = TRUE, recursive = FALSE)
  
Candida_instruments <- fread(file="/Volumes/Passport/119Bacs_addedP/addedPgenus.CandidatusSoleaferrea.id.11350.summary.txt.csv",header = T, sep = ",")

for (k in files) {
  TOCOLOC <- read.table(k, header = T, sep="\t")
  nameOfTest <- str_replace(str_split(k, "/", n = Inf, simplify = TRUE)[[8]], ".tab","")
  
  #Just need a subset of file that contains SNPs with significance of 5E-N for both disease and phenotype
  SNPlist <- subset(TOCOLOC, select=c(topSNP, topSNP_bp, topSNP_chr)) 
  
  #For each SNP, search through phenotype summary stats to get 1MB of variants b4 and after SNP, this is called  region
  for(i in 1:length(SNPlist$topSNP)) {
    thisSNP <- SNPlist$topSNP[i]
    thisChr <- SNPlist$topSNP_chr[i]
    thisBP <- SNPlist$topSNP_bp[i]
    thisBpLow <- thisBP - 1000000
    thisBpHigh <- thisBP + 1000000
    getRegions <- subset(Candida_instruments, chr == thisChr & bp > thisBpLow & bp < thisBpHigh)
    print(paste("Count is:", i))
    print(paste("Size of region for ",thisSNP," in ",thisChr,"_",thisBP, " : ", nrow(getRegions), sep=""))
    dir.create(paste("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/coloc_mrQTL_regions/",nameOfTest,sep=""))
    write.table(getRegions, paste("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/coloc_mrQTL_regions/",nameOfTest,"/",thisSNP,"_region.tab",sep = ""), quote = F, sep = "\t", row.names = F)
  }
}


#(PD vs Candida)
##Bayesian colocalisation analysis

############Excluding some columns and changing name of columns of region files############
library("dplyr")
masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/regionsV2"
files <- list.files(path=masterDir, pattern="*.tab", full.names = TRUE, recursive = FALSE)

for (i in files) {
  #/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/regions/rs489567_region.tab
  data <- fread(file=i,header=T,sep="\t")
  df <- select(data, `chr:pos`, rsID, eff.allele, ref.allele, `beta.x`, SE, N, Ncohorts, pDerived)
  names(df) <- c("SNP","rsID", "A1", "A2", "beta", "SE", "N", "Ncohorts","P-value")
  write.table(df, i, quote = F, sep = "\t", row.names = F)
}
##############################################################################################################

############reading one type of Dx sum stats and renaming col so that merging is possible############

dataPD <- fread(file="/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/META5_all_with_rsid.txt",header=T,sep="\t")
names(dataPD)[17] <- "rsID"
write.table(dataPD, "/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/META5_all_with_rsid.tab", quote = F, sep = "\t", row.names = F)

##############################################################################################################

############reading one type of Dx sum stats but no change of col bc rsID column name already exists for merging############
nalls_PD <-fread("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/nallsEtAl2019_excluding23andMe_allVariants.tab",header=T,sep="\t")
##############################################################################################################

nalls_PD$SNP <- gsub("chr","",nalls_PD$SNP) #substituing a pattern of characters in a specific column in Dx sumstats

############Merging all info on candida data  with even more info from the Dx sum stats by rsID############

masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/regionsFullStats"
files <- list.files(path=masterDir, pattern="*.tab", full.names = TRUE, recursive = FALSE)

for (i in files) {
  print(paste("Working on...",i))
  data <- fread(file=i,header=T,sep="\t")
  print("here")
  df <- merge(data, dataPD, by="rsID")
  print("here")
  write.table(df, i, quote = F, sep = "\t", row.names = F)
}


##############################################################################################################

############Running the colocalization############

#We subset data frames and rename them and let the function coloc.abf do all the work
SNPlist <- fread(file = paste("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/hitsBelow5E-5ForPD/hitsBelow5E-4forCandida/Below5E-4.varsOf_PD_and_CandidatusSol_sumstats.tab"), header = T, sep = "\t")

regions <- c(SNPlist$ID)

for (i in regions) {
  #/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/regions
  thisRegion <- i
  print(i)
  data <- fread(file = paste("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/regionsFullStats/",thisRegion, "_region.tab", sep = ""), header = T, sep = "\t")
  #print(names(data))
  df2_Candida_coloc <- data[,c("SNP","beta","P-value.x")]
  df1_PD_coloc <- data[,c("SNP","b","p","Freq1")]
  df2_Candida_coloc$var <- var(df2_Candida_coloc$beta)
  df1_PD_coloc$var <- var(df1_PD_coloc$b)
  names(df1_PD_coloc) <- c("snp","beta", "pvalues","MAF","varbeta")
  names(df2_Candida_coloc) <- c("snp","beta", "pvalues","varbeta")
  
  #Remove duplicates in dataset 1
  df1_PD_coloc <- df1_PD_coloc[!duplicated(df1_PD_coloc$snp), ]
  df2_Candida_coloc <- df2_Candida_coloc[!duplicated(df2_Candida_coloc$snp), ]
  
  df1_PD_coloc$type <- "cc"
  df2_Candida_coloc$type <- "quant" 
  df1_PD_coloc$s <- 0.0590704464
  df2_Candida_coloc$sdY <- 1
  results <- coloc.abf(df1_PD_coloc, df2_Candida_coloc, MAF = NULL, p1 = 1e-02, p2 = 1e-02, p12 = 1e-02)
  #results <- coloc.abf(df1_PD_coloc, df2_Candida_coloc, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05) 
  sink(file = paste("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/PD_coloc/",thisRegion,"_resultsSummary.txt",sep = ""))
  print(results$summary)
  sink()
  write.table(results$results, file = paste("/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/PD_coloc/",thisRegion,"_resultsSnps.txt",sep = ""), quote = F, sep = ",", row.names = F)
}
##########################################################################################################################
#####POST COLOCALIZATION: FINDING COLOCALIZED REGIONS THAT EXPLAIN PHENOTYPE AND DX ETIOLOGIES VIA SIMILAR MECHANISMS#####
##########################################################################################################################

masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/Colocolization/PD_coloc"
files <- list.files(path=masterDir, pattern="*Summary.txt", full.names = TRUE, recursive = FALSE)

for (i in files) {
  #print(paste("Working on...",i))
  data <- fread(file=i,header=T,sep=" ")
  #print(paste("The PP.H4.abf for this region is:",data$`PP.H4.abf`))
  print(as.double(data$PP.H4.abf))
  if (as.double(data$PP.H0.abf) > 0.95) {
    print("PP.H0.abf HIT")
  }
  if (as.double(data$PP.H1.abf) > 0.95) {
    print("PP.H1.abf HIT")
  }
  if (as.double(data$PP.H2.abf) > 0.95) {
    print("PP.H2.abf HIT")
  }
  if (as.double(data$PP.H3.abf) > 0.95) {
    print("PP.H3.abf HIT")
  }
  if (as.double(data$PP.H4.abf) > 0.95) {
    print("PP.H4.abf HIT")
  }
}

##########################################################################################################################





