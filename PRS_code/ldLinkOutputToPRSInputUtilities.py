#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 22:28:55 2021

@author: rodrigosandon
"""
import pandas as pd
import csv


class PostLDLinkUtilities:
    
    def __init__(self, file_path):
        self.file_path = self
        
        
csv = "/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned/exampleOfLDLinkOutput.csv"

df = pd.read_csv(csv)

lstOfKeptVariants = list(df.loc[df["Details"] == "Variant kept."]["RS Number"])#[specifying what rows][get which column?]

csvOfFileBeingPruned = "/Volumes/Passport/MiBioGen_filtered/addedPgenus..Clostridiuminnocuumgroup.id.14397.summary.txt.csv"

df2 = pd.read_csv(csvOfFileBeingPruned)

lstOfKeptVariantsGoodFormat = list(df2.loc[df2["rsID"].isin(lstOfKeptVariants)]["chr"].astype(str) + ":" + df2["bp"].astype(str))
lstOfKeptEffectAlleles = list(df2.loc[df2["rsID"].isin(lstOfKeptVariants)]["eff.allele"])
lstOfKeptBetaValues = list(df2.loc[df2["rsID"].isin(lstOfKeptVariants)]["beta"])

cleanedList1 = [i for i in lstOfKeptVariantsGoodFormat if str(i) != 'nan']
cleanedList2 = [i for i in lstOfKeptEffectAlleles if str(i) != 'nan']
cleanedList3 = [i for i in lstOfKeptBetaValues if str(i) != 'nan']

data = {'rsID' : cleanedList1, 'eff.allele' : cleanedList2, 'beta' : cleanedList3}

newDf = pd.DataFrame(data, index=None, columns=None)

###Seeing if there are duplicates
len(cleanedList1) != len(set(cleanedList1)) #<--- to see if there duplicates in the file, if so, build code to prune it

#actually converts it to text file
newDf.to_csv('/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned/exampleScoreFile.txt', header=None, index=None, sep=' ', mode='a')



