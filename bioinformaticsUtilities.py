#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 07:15:28 2021

@author: rodrigosandon
"""

import pandas as pd
import time

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
                               delim_whitespace=True)
        
        readFile.to_csv(outputNamePath, index=None)
        
    def readHead (csvPath):
        df = pd.read_csv(csvPath)
        print(df.head())
        print(df.size)
        return df.head()
    
    def matchAndFill (output, toFill) :
        df1 = pd.read_csv(output)
        df2 = pd.read_csv(toFill)
        
        
        
    
#genomicUtilities.txtFileToCSVFile("/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/nallsEtAl2019_excluding23andMe_allVariants_singleSpaced.txt", "/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/nalls_dataset.csv")
#csvPath = "/Volumes/Passport/NIHSummer21/G_Lactobacillus_HB_allchr.txt"
#df = genomicUtilities.readHead(csvPath)
#genomicUtilities.txtFileToCSVFile("/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/hg19_positions.txt", "/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/hg19_positions.csv")
genomicUtilities.matchAndFill("/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/hg19_positions.csv", "/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/allMB_variants.csv")

end = time.time()
print(end - start)
