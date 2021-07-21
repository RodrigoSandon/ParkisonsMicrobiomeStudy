#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 07:15:28 2021

@author: rodrigosandon
"""

import pandas as pd
import time
import csv
import os

start = time.time()

class genomicUtilities:
    
    def __init__(self, file_path):
        self.file_path = self
    
    def txtFileToCSVFile (txtFilePath, outputNamePath):
        txtin = open(txtFilePath, "rt")
        csvout = open(outputNamePath, "wt")

        for line in txtin:
            csvout.write(' '.join(line.split()))
            csvout.write(' \n')

        txtin.close()
        csvout.close()
        
        
        readFile = pd.read_csv(txtFilePath, 
                               delim_whitespace=True,  error_bad_lines=False)
        
        readFile.to_csv(outputNamePath, index=None)
        
        
    def readHead (csvPath):
        df = pd.read_csv(csvPath)
        print(df.head())
        print(df.size)
        return df.head()
        
    def makeTxtFiles(path):
        for root, dirs, files in os.walk(path):
            #print(dirs)
            newFile = os.path.join(root, "prunedVariants.txt")
            file = open(newFile, "w")
            file.close()
    
    def processMultipleFolders (path):
        count = 0
        for root, dirs, files in os.walk(path):
            for name in files:
                curr = os.path.join(root, name)
                if name.endswith("prunedVariants.txt"):
                    print(name)
                    count = count + 1
                    splt = curr.split(".")
                    toReplace = splt[0] + str(".csv")
                    #print(toReplace)
                    genomicUtilities.txtFileToCSVFile(curr, toReplace)
        print("Files processed: ", count)
        
    def checkMatchesInTwoCSVs(csv1Path, csv2Path): #returns list of rsIDs that are found in both CSVs
        arrOfMatches = []
        numberOfMatches = 0
        
        csv1 = pd.read_csv(csv1Path)
        csv2 = pd.read_csv(csv2Path)
        
        listOfCSV1 = csv1["rsID"].tolist() #hardcoded column
        #print(listOfCSV1)
        #print(len(listOfCSV1))
        listOfCSV2 = csv2["rsID"].tolist()
        #print(listOfCSV2)
        #print(len(listOfCSV2))
        
        for i in listOfCSV1:
            if i in listOfCSV2:
                arrOfMatches.append(i)
                numberOfMatches = numberOfMatches + 1
        
        #print(numberOfMatches)
        #print(arrOfMatches)
        
        return arrOfMatches
    
    def outputCSVFromMatches (filetoExtractFrom, csv2Path, outCSV):
        matches = genomicUtilities.checkMatchesInTwoCSVs(filetoExtractFrom, csv2Path)
        #print(matches)
        #print(len(matches))
        df = pd.read_csv(filetoExtractFrom)
        
        subsetDf = []
        
        for row in range(len(df.index)):
            for col in range(len(df.columns)):
                #print(row,"",col)
                if df.iloc[row][col] in matches:
                    #print(df.iloc[row][col])
                    #print(df.loc[col].tolist()) <--- works, each series to list, one column is an element
                    #print(df.loc[col]) gets row
                    #print(type(df.loc[col])) <---- series
                    subsetDf.append(df.iloc[row].tolist()) #this is how you append to row as list to a list
                    #print(df.iloc[row].tolist())
        #print(subsetDf)
        new_df = pd.DataFrame(subsetDf, columns = df.columns.tolist())
        new_df.to_csv(outCSV, index =False, header = True)
        
    def CSVtoTXT (csv_file, txtFile):

        with open(txtFile, "w") as my_output_file:
            with open(csv_file, "r") as my_input_file:
                [my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
            my_output_file.close()
            
    def extractOnlyCertainCols (csvPath, csvOutPath):
        print("easier if manually")
        
    def dropDuplicates (csvPath, pathOut):
        
        df = pd.read_csv(csvPath)
        
        dupRows = df[df.duplicated(subset = ["chrID"])]
        #print(dupRows)
        print(len(dupRows))
        
        noDups_df = df.drop_duplicates(subset = ["chrID"])
        
        noDups_df.to_csv(pathOut, header = False, index = False)

###Functions for supplementary informartion
    def CSVofFileNamesInFolder (pathToFolder, pathOfCSVOutFile):
        bacteriaNames = []
        numberOfSamplesInEachGenus = []
        count = 0
        for root,dirs,files in os.walk(pathToFolder):
            for name in files:
                if name.endswith(".csv"):
                    splitt = name.split(".")
                    try:
                        nameOfBacteria = splitt[1]
                        if nameOfBacteria == '':
                            nameOfBacteria = splitt[2]
                    except IndexError:
                        pass
                    #print(nameOfBacteria)
                    print(count)
                    df1 = pd.read_csv(os.path.join(root,name))
                    numberOfSamplesInEachGenus.append(len(df1))
                    if len(nameOfBacteria) > 0 and nameOfBacteria != "_genus": #excluding some weird names from count
                        bacteriaNames.append(nameOfBacteria)
                        count = count + 1
                
        #print(bacteriaNames)
        print(count)
        data = {'Genus' : bacteriaNames, 'Samples Avaliable' : numberOfSamplesInEachGenus}
        df = pd.DataFrame(data)
        df.to_csv(pathOfCSVOutFile, header = True, index = False)
        """with open(pathOfCSVOutFile, 'w') as myfileToWriteOn:
            wr = csv.writer(myfileToWriteOn, quoting=csv.QUOTE_ALL)
            try:
                wr.writerow(bacteriaNames)
            except TypeError:
                print("TypeError!")"""
                

genomicUtilities.CSVofFileNamesInFolder("/Volumes/Passport/MiBioGen_filtered", "/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen/genusSuppInfo.csv")

###############################################################
#For creating csv files in multiple folders

"""path = "/Volumes/Passport/MiBioGen_QmbQTL_summary_genus/"

toCSVSoFar = []
for root, dirs, files in os.walk(path):
    for name in files:
        curr = os.path.join(root, name)
        if curr.endswith(".csv"):
            toCSVSoFar.append(curr)

for root, dirs, files in os.walk(path):
    for name in files:
        curr = os.path.join(root, name)
        if curr.endswith(".txt") and curr + str(".csv") not in toCSVSoFar:
            #splt = curr.split(".")
            #toReplace = splt[0] + str(".csv")
            toReplace = curr + str(".csv")
            print(toReplace)
            genomicUtilities.txtFileToCSVFile(curr, toReplace) #<--- commented out once done"""


###############################################################
#genomicUtilities.matchAndFill("/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/hg19_positions.csv", "/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/allMB_variants.csv")

end = time.time()
print(end - start)
