install.packages("LDlinkR")
library("LDlinkR")
library("hash")
library("stringr")
library(dplyr)
library(biomaRt)
library(tidyverse)
library(sjmisc)
csv = "/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/ld_linkAccInput.csv"
#output is the 2 files from the api into correct folders
#mb_markers1 <- mbMarkers_martins %>% select(4:5)
#LDmatrix(snps, pop = "CEU", r2d = "r2", token = NULL, file = FALSE)
ld_input = read.csv(csv)

mbDict <- hash()

for (col in 1:ncol(ld_input)) {
  colName <- paste("chr",col,sep="")
  #print(colName)
  lst <- list()
  for (row in 1:49) {
    print(colName)
    print(row)
    element <- ld_input$colName[row]
    print(element)
    if (str_contains(ld_input$colName[row], '.') == FALSE) {
      lst <- c(lst, ld_input$colName[row])
    }
  }
  #lst <- ld_input %>% select(col) #select from a give data.frame, all the elements in the column within select()
  mbDict[[colName]] <- lst
}

#problem -> trouble creating a dictionary in R so that that I have a list of SNPS per chromosome

for (v in values(mbDict["chr8"])) {
  print(str_contains(v, "."))
  if (str_contains(v, ".") == TRUE) {
    print(v)
  }
}
