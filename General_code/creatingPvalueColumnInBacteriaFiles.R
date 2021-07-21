#### Do calculation of P-value first to all bacterias

addPValColumnAndFilterP <- function(csvToRead, newDirToPlaceNewColumnCSV) { # make sure to ad "/" at the end of arg2
  start.time <- Sys.time()
  temp <- read.csv(csvToRead, header = TRUE, sep = ",") 
  temp$pDerived <- 2*pnorm(-abs(temp$beta/temp$SE))
  print("Here")
  temp1 <- temp[!(temp$pDerived > 0.05),] #omitting p values over 0.05
  print("Here")
  splitt <- strsplit(csvToRead, "/")
  #print(splitt)
  newNameForDf <- paste("addedP", splitt[[1]][5], sep = "") #5 is where the name of file is found in list
  #print(newNameForDf)
  newPath <- paste(newDirToPlaceNewColumnCSV,newNameForDf,sep = "")
  write.csv(temp1, newPath, row.names = FALSE)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

masterDir <- "/Volumes/Passport/MiBioGen_QmbQTL_summary_genus"
newDirToPlaceNewColumnCSV = "/Volumes/Passport/MiBioGen_filtered/"
files <- list.files(path=masterDir, pattern="*.csv", full.names = TRUE, recursive = FALSE)
for (i in files) {
  print(i)
  addPValColumnAndFilterP(i, newDirToPlaceNewColumnCSV)
}
