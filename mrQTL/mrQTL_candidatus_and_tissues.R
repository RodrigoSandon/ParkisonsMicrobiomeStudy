
library(data.table)
library(stringr)

masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/MR_QTL"
files <- list.files(path=masterDir, pattern="*.msmr", full.names = TRUE, recursive = FALSE)

listOfTissueNumberVariants <- c()

for (i in files) {
  print(paste("Working on:", i))
  Candida <- fread("/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/CandidatusSoleaferrea.id.11350.toscore.txt", header =F)
  colnames(Candida) <- c("chr_pos", "A1", "beta")
  temp <- fread(i, header =T)
  listOfTissueNumberVariants <- c(listOfTissueNumberVariants, nrow(temp))
  temp$chr_pos <- paste(temp$topSNP_chr, temp$topSNP_bp, sep = ":") #make temp data.frame of merged columns from tissue sample
  total <- merge(Candida, temp, by="chr_pos")
  #newName = str_replace(i,"msmr","txt")
  #write.table(total, file = paste(i, ".txt", sep=""), quote = F, row.names = F, sep = "\t")
}


masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/results"
files2 <- list.files(path=masterDir, pattern="*.txt", full.names = TRUE, recursive = FALSE)
#/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/results/... 7
for (i in files2) {
  print(paste("Working on:", i))
  result <- fread(i, header =T)
  fileName <- str_split(i, "/", n = Inf, simplify = TRUE)[[8]]
  #print(fileName)
  newFileName <- str_replace(paste("/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/results_csv/",fileName,sep=""),"msmr.txt","csv")
  print(paste("Find in:", newFileName))
  write.csv(result, file = newFileName, row.names = FALSE)
}

masterDir <- "/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/results"
files2 <- list.files(path=masterDir, pattern="*.txt", full.names = TRUE, recursive = FALSE)
#/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/results/... 7

count <- 0
for (i in files2) {
  count <- count + 1 #ahead bc indexing a vector starts w 1
  print(paste("Working on:", i))
  result <- fread(i, header =T)
  hits <- result[result$p_SMR_multi < 0.05/listOfTissueNumberVariants[count], ]
  print(0.05/listOfTissueNumberVariants[count])
  print(nrow(hits))
  fileName <- str_split(i, "/", n = Inf, simplify = TRUE)[[8]]
  newFileName <- str_replace(paste("/Volumes/T7Touch/NIHSummer2021/Code/MR_QTL_analysis/p_SMR_multi_belowBonferroni/",fileName,sep=""),"msmr.txt","csv")
  write.csv(hits, file = newFileName, row.names = FALSE)
}
