#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 14:48:17 2021

@author: rodrigosandon
"""
import ldLinkPipeline
import time
import os

parent_dir = "/Volumes/Passport/MiBioGen_filtered"

dir_for_results = "/Volumes/Passport/MiBioGen_LDPruning"

listOfCSVPathsToProcess = []

for root,dirs,files in os.walk(parent_dir):
    for name in files:
        listOfCSVPathsToProcess.append(os.path.join(root, name))

###THE PIPELINE
for i in reversed(listOfCSVPathsToProcess): #i is the og bac sum stats file

    print("Processing %s...."%(i))
    start = time.time()
    
    df = ldLinkPipeline.ldLinkUtilities.CSVtoDataFrame(i)
    
    ldLinkPipeline.ldLinkUtilities.makeChrIDColumn(df)
    
    #specified column to a list
    listOfChrIds = list(df["chrID"])
    #prepare the dictionary
    dOfChrIds = {}
    
    ldLinkPipeline.ldLinkUtilities.chrIDColumnToDictionary(listOfChrIds,dOfChrIds)
    
    listOfStringUnixInputs = ldLinkPipeline.ldLinkUtilities.listOfListsToListOfStrings(ldLinkPipeline.ldLinkUtilities.dictToList(dOfChrIds))
    
    dir_toPlaceLDLinkOutput = "/Volumes/Passport/MiBioGen_LDPruning"
    
    folderNameToPutResultsIn = ldLinkPipeline.ldLinkUtilities.makeAndGetFolderNameForBac(i, dir_toPlaceLDLinkOutput)
    
    ldLinkPipeline.ldLinkUtilities.runLDLinkCommandForAListOfInputs(listOfStringUnixInputs, dir_for_results, folderNameToPutResultsIn)
    
    end = time.time()
    print("Time taken to complete bacteria file was  %s seconds"%(end - start))