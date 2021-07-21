
addPColumn <- function(csvToRead, newDirToPlaceNewColumnCSV) { # make sure to ad "/" at the end of arg2
  start.time <- Sys.time()
  temp <- read.csv(csvToRead, header = TRUE, sep = ",") 
  temp$pDerived <- 2*pnorm(-abs(temp$beta/temp$SE))
  print("Here")
  splitt <- strsplit(csvToRead, "/")
  newNameForDf <- paste("addedP", splitt[[1]][5], sep = "") #5 is where the name of file is found in list
  newPath <- paste(newDirToPlaceNewColumnCSV, newNameForDf,sep = "")
  write.csv(temp, newPath, row.names = FALSE)
  end.time <- Sys.time()
  print(end.time - start.time)
}

masterDir <- "/Volumes/Passport/MiBioGen_QmbQTL_summary_genus"
newDirToPlaceNewColumnCSV = "/Volumes/Passport/119Bacs_addedP/"
files <- list.files(path=masterDir, pattern="*.csv", full.names = TRUE, recursive = FALSE)

### Loop through all csv file paths
for (i in files) {
  print(paste("Working on:", i))
  if (grepl("unknown",i, fixed = TRUE) == FALSE) { #meaning the string doesn't contain unknow, then process
    addPColumn(i, newDirToPlaceNewColumnCSV)
  }
}


#csvToRead = "/Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/nalls_onlyRsIDs.csv"
#temp <- fread(csvToRead, header = T)
#temp$p <- 2*pnorm(-abs(temp$beta/temp$SE))
#write.table(temp, file = "/Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/nalls_onlyRsIDs.csv", sep = ",",
#            row.names = FALSE, quote = FALSE)