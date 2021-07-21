library("data.table")
library("dplyr")
library("tidyr")
library("plyr")
#install.packages("rjson")
library("rjson")

#CSV to text PD sum stats
PD_sumstats = "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/nalls_dataset_formatted.csv"
PD <- fread(PD_sumstats, header = T)

write.table(PD, file = "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/nalls_dataset_formatted.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

#CSV to text BAC sum stats

MB_sumstats = "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/bacEx_datasetFormatted.csv"
MB <- fread(MB_sumstats, header = T)

write.table(MB, file = "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/bacEx_datasetFormatted.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

#Need to use unpruned ones, I forgot

MB_sumstats = "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/bacExample_dataset_formatted.csv"
MB <- fread(MB_sumstats, header = T)

write.table(MB, file = "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/bacExample_dataset_formatted.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

#If you want to compute the genetic correlation between schizophrenia and bipolar disorder, type these commands
# Download Data

