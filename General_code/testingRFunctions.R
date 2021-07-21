
library(data.table)
bim <- fread("/Volumes/T7Touch/NIHSummer2021/Data/PLINK_analysis/IMPUTED.SOFTCALLS.Demo.bim", header =F)
colnames(bim) <- c("column1", "column2", "column3", "column4", "column5", "column6")
bim$columngood <- paste(bim$column1, bim$column4, sep = ":")
outPut <- bim[,c("column1", "columngood", "column3", "column4", "column5", "column6")]
#outPut2 <- unique(outPut, by ="columngood")

duplicated(outPut$columngood)
outPut$columngood[duplicated(outPut$columngood)]
duplicates <- outPut$columngood[duplicated(outPut$columngood)]
write.table(outPut, file = "/Volumes/T7Touch/NIHSummer2021/Data/PLINK_analysis/temp.bim", quote = F, row.names = F, col.names = FALSE, sep = "\t") 

df <- as.data.frame(outPut)

#dt1 <- data.table(column1 = c(), columngood = c(), column3 = c(), column4 = c(), column5 = c(), column6 = c())
newdf <- c()
count <- 0
for (i in df$columngood) { #df and outPut are the same table
  #print(paste(i,"out of",13085227), sep = " ")
  if (i %in% duplicates) {
    print(i)
    print(list(outPut[count,]))
    newdf <- c(newdf, list(outPut[count,]))
  }
  count = count + 1
}

hopedf <- data.frame(matrix(unlist(newdf), nrow=length(newdf), byrow=TRUE),stringsAsFactors=FALSE)


### okay so i can treat an R data.frame as a matrix

######RUN HERE: dropping rows that contain ids that are not in the duplicates hash table (h) from df copy dataframe (df_only_dups) 
df_only_dups <- as.data.frame(df)

for (i in 1:nrow(df)){
    print(i)
    splitt <- str_split(df[i,2],":", simplify = TRUE)
    chrID <- splitt[1,1]
    if (chrID != '3' || chrID != '5') { #check if the id is in 3 or 5, if so, don't process
      chrIDint <- as.integer(chrID)
      if (!(df[i,2] %in% values(h,keys=chrIDint))) { # if not in the list of the key chrIDint
        print(paste("dropping....",df[i,2], sep = " "))
        df_only_dups <- df_only_dups[-c(i),] #then drop entire row
      }
    }
}
######STOP

##### more efficient way^
#install.packages("stringr")
library(stringr)
library(hash)

h <- hash()

for (i in duplicates) {
  splitList <- str_split(i,":", simplify = TRUE)
  key <- splitList[1,1]
  value <- i
  
  if (!(has.key(key, h))) { #if no key value pair existing
    h[key] <- list()
    h[key] <- append(h[[key]], value)
  } else {
    h[key] <- append(h[[key]], value) #if exist already
  }
}
#first turn list of duplicated ids into their respective chromosome #s
#now have has

#counting all values in hash
countOfids_inH <- 0
chrs <- c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,4,6,7,8,9) #chrs 3 and 5 dont have duplicates
for (i in chrs) {
  lenOfKey <- length(values(h,keys=i))
  countOfids_inH <- countOfids_inH + lenOfKey
}


#way to slow, took like 6hrs to get to 1350, at this rate --> 