#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:11:13 2021

@author: rodrigosandon

"""

import pandas as pd
import time
import os

class ldLinkUtilities:
    
    def __init__(self, file_path):
        self.file_path = self
        
    def appendTo2DArrayOrMakeNewArrayIn2DArray (TwoDarray, valueToAppend):
        added = False
        for i in range(len(TwoDarray)):
            if len(TwoDarray[i]) < 1000:
                TwoDarray[i].append(valueToAppend)
                added = True
        if added == False: #meaning array or all arrays were full, make new array and add onto there
            TwoDarray.append([]) #should append to last of list?
            TwoDarray[-1].append(valueToAppend)
            added = True
            
    def hasKey (d, key):
        if key in d:
            return True
        else:
            return False
        
    def numberOfElementsInDict (d):
        count = 0
        for key in d:
            for lst in d[key]: #d[key] is the entire list, within this list is x amount of lists
                count += len(lst)
        #print(count)
                
    def listOfListsToListOfStrings (listToChangeElementsOf):
        
        for i in range(len(listToChangeElementsOf)):
            tempString = ""
            tempString += "\""
            for j in range(len(listToChangeElementsOf[i])): #now at a specific list within the list
                if isinstance(listToChangeElementsOf[i][j], float) == False:
                    if j < len(listToChangeElementsOf[i]) - 1: # if not at the end
                        #print(listToChangeElementsOf[i][j])
                        tempString += "".join((listToChangeElementsOf[i][j],"\\n")) #separator is "\n"
                    else: #if at the end
                        tempString += listToChangeElementsOf[i][j]
            tempString += "\""
            listToChangeElementsOf[i] = tempString #actually replacing the elements for a string
        return listToChangeElementsOf
    
    ###DEPRICATED#####
    def dictToString (d, sep):
        
        result = ""
        for key in d:
            for string in d[key]:
                if string != d[key][-1]: #if not last element
                    result += string + sep #sep should be " " for a unix command
                else:
                    result += string
        return result
    
    def TwoDimArrToString (lst, sep):
        
        result = ""
        for string in lst:
            if string != lst[-1]: #if not last element
                result += string + sep #sep should be " " for a unix command
            else:
                result += string
        return result
    
    ######VERY SLOW########
    def replaceChridForRsid (newD, dfasReference):
        for key in newD:
            start = time.time()
            for lstIdx in range(len(newD[key])): #now at indv chrIDs
                for chrIdIdx in range(len(newD[key][lstIdx])):
                    chrId = newD[key][lstIdx][chrIdIdx]
                    #print(chrId)
                    rsId = dfasReference.loc[dfasReference["chrID"] == chrId]["rsID"] ####CRUCIAL PART-->finds rsID equivalent of chrID
                    #print(rsId.item())
                    newD[key][lstIdx][chrIdIdx] = rsId.item() #the count1 and count2 serve as indices
                    #print(newD[key][lstIdx][chrIdIdx]) 
            end = time.time()
            time_taken = end - start
            print(f'Done replacing variants in {key}, time taken : {time_taken}.')
        return newD
    
    ####NEW ONE FOR ^#######
    def replaceChridIn2DList (lst, dfasReference):
        for row in range(len(lst)): #now at indv chrIDs
            start = time.time()
            for col in range(len(lst[row])):
                chrId = lst[row][col]
                rsId = dfasReference.loc[dfasReference["chrID"] == chrId]["rsID"] ####CRUCIAL PART-->finds rsID equivalent of chrID
                lst[row][col] = rsId.item()
            end = time.time()
            time_taken = end - start
            print(f'Done replacing variants in row {row}, time taken : {time_taken}.')
        return lst
    
    def dictToList (d):
        masterlist = []
        for key,value in d.items():
            for lsts in value:
                masterlist.append(lsts)
        return masterlist
    
    def CSVtoDataFrame (csv):
        df = pd.read_csv(csv)
        return df
    
    def makeChrIDColumn (df):
        df.insert(3,"chrID", "chr" + df["chr"].apply(str) + ":" + df["bp"].apply(str), allow_duplicates=True) 
        
    def chrIDColumnToDictionary (lst,dOfChrIds): #Data structure of dict: {key : [[],[],...[]], key: [[],[],...[]], ... key: [[],[],...[]]}
        for i in range(len(lst)):
            keyName = lst[i].split(":")[0]
            
            if ldLinkUtilities.hasKey(dOfChrIds, keyName) == True: #most likely it will already exist
                ldLinkUtilities.appendTo2DArrayOrMakeNewArrayIn2DArray (dOfChrIds[keyName], lst[i])
            else: #will only happen the first time for each chr, should never be over 1000 bc just started
                dOfChrIds[keyName] = [[]]
                dOfChrIds[keyName][0].append(lst[i])
    
    def makeAndGetFolderNameForBac (ogCSVPathOfBacVariants, parent_dir):
        nameOfBacteria = ""
        splitt = ogCSVPathOfBacVariants.split(".")
        try:
            nameOfBacteria = splitt[1]
            if nameOfBacteria == '':
                nameOfBacteria = splitt[2]
        except IndexError:
            pass
        
        print(nameOfBacteria)
        folderPath = os.path.join(parent_dir, nameOfBacteria)
        try:
            os.mkdir(folderPath)
        except OSError:
            pass
        return nameOfBacteria
        
    def runLDLinkCommandForAListOfInputs(listOfStringUnixInputs, root_dir, folderNameToPutResultsIn):
        count = 1
        
        ##NEW FEATURE: Identify where it left off (at what element in the list of string inputs)
        listOfCurrentOutputFilesInBacteriaFile = []
        #print(os.path.join(root_dir, folderNameToPutResultsIn))
        for root, dirs, files in os.walk(os.path.join(root_dir, folderNameToPutResultsIn)):
            for name in files:
                listOfCurrentOutputFilesInBacteriaFile.append(name)
        
        for i in listOfStringUnixInputs: # the list is the string calls
            try:
                chrom = i.split("\"")[1].split("\\n")[0].split(":")[0] #both of these try and catach do same thing, just that one (first one) is for the exception when the list of stringinput includes only one  SNP
            except KeyError:
                selectEleToGetChrom = i.split("\\n")[1]
                chrom = selectEleToGetChrom.split(":")[0]
            #selectEleToGetChrom = i.split("\\n")[1]
            #chrom = selectEleToGetChrom.split(":")[0]
            filename = "%s_%s.csv"%(chrom,count)
            print(filename)
            if filename in listOfCurrentOutputFilesInBacteriaFile:
                count += 1
                continue
            else:
                #other parameters that can be changed: pop, r2, maf
                cmd = "curl -k -H \"Content-Type: application/json\" -X POST -d '{\"snps\": %s, \"pop\": \"EUR\", \"r2_threshold\": \"0.3\", \"maf_threshold\": \"0.01\"}' 'https://ldlink.nci.nih.gov/LDlinkRest/snpclip?token=4dd19de0dba7' >> %s/%s/%s"%(i, root_dir, folderNameToPutResultsIn, filename)
                count += 1
                
                os.system(cmd)
                
###Using another token^: joe.gulags@gmail.com : 4dd19de0dba7
################################################################################################################

#Enter entire MiBioGen_filtered folder that contains all 119 bacterial sum stats

"""parent_dir = "/Volumes/Passport/MiBioGen_filtered"

listOfCSVPathsToProcess = []

for root,dirs,files in os.walk(parent_dir):
    for name in files:
        listOfCSVPathsToProcess.append(os.path.join(root, name))

###THE PIPELINE
for i in listOfCSVPathsToProcess: #i is the og bac sum stats file
    print("Processing $s...."%(i))
    start = time.time()
    
    df = ldLinkUtilities.CSVtoDataFrame(i)
    
    ldLinkUtilities.makeChrIDColumn(df)
    
    #specified column to a list
    listOfChrIds = list(df["chrID"])
    #prepare the dictionary
    dOfChrIds = {}
    
    ldLinkUtilities.chrIDColumnToDictionary(listOfChrIds,dOfChrIds)
    
    listOfStringUnixInputs = ldLinkUtilities.listOfListsToListOfStrings(ldLinkUtilities.dictToList(dOfChrIds))
    
    dir_toPlaceLDLinkOutput = "/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned/LDLink_Output"
    
    folderNameToPutResultsIn = ldLinkUtilities.makeAndGetFolderNameForBac(i, dir_toPlaceLDLinkOutput)
    
    ldLinkUtilities.runLDLinkCommandForAListOfInputs(listOfStringUnixInputs, folderNameToPutResultsIn)
    
    end = time.time()
    print("Time taken to complete bacteria file was  %s seconds"%(end - start))

################################################################################################################


################################################################################################################
#load csv
csv1 = "/Volumes/Passport/MiBioGen_filtered/addedPgenus..Clostridiuminnocuumgroup.id.14397.summary.txt.csv"
df = pd.read_csv(csv1)

#create a new column that is the chrID (format required for LDLink)
df.insert(3,"chrID", "chr" + df["chr"].apply(str) + ":" + df["bp"].apply(str), allow_duplicates=True) 
#^this is for chrIDs, but really want rsIDs

#make dictionary of chrIDs, if over 1000, make a new key-value pair for that chr
listOfChrIds = list(df["chrID"])
dOfChrIds = {}

###from the list of the column df, arrange this into a dict
for i in range(len(listOfChrIds)):
    keyName = listOfChrIds[i].split(":")[0]
    
    if hasKey(dOfChrIds, keyName) == True: #most likely it will already exist
        appendTo2DArrayOrMakeNewArrayIn2DArray (dOfChrIds[keyName], listOfChrIds[i])
    else: #will only happen the first time for each chr, should never be over 1000 bc just started
        dOfChrIds[keyName] = [[]]
        dOfChrIds[keyName][0].append(listOfChrIds[i])
        
#At this point, the way its organized with a dict doesn't matter , so the output from this pyfile can be 1 huge list of strings now
#    instead of a dictionary, in the shell script, ill do the for loop 

#convert dict to 2d list
TwoDListOfAllChrs = dictToList(dOfChrIds)

print("Here")
#replace all chrIDs for rsIDs
#newDictOfrsIDs = replaceChridForRsid (dOfChrIds, df)
#rsIDs2dList = replaceChridIn2DList (TwoDListOfAllChrs, df) <----- NO LONGER NEEDED TO REPLACE CHRID W/ RSID

#converting list of lists to list of strings
#listOfStringUnixInputs = convertListsInListIntoString(rsIDs2dList[20:]) 
listOfStringUnixInputs = listOfListsToListOfStrings(TwoDListOfAllChrs)  #the convertListsInListIntoString affects both input and output
print("Here")

#####this is getting the ldlink output files#####this is getting the ldlink output files#####this is getting the ldlink output files

y = ["\"chr18:6344587\\nchr18:6343719\\nchr18:6343212\\nchr18:6348227\\nchr18:6354297\\nchr18:6343636\\nchr18:6345363\\nchr18:6345071\\nchr18:6348470\\nchr18:6348553\\nchr18:6348028\""
     , "\"chr18:6344587\\nchr18:6343719\\nchr18:6343212\\nchr18:6348227\\nchr18:6354297\\nchr18:6343636\\nchr18:6345363\\nchr18:6345071\\nchr18:6348470\\nchr18:6348553\\nchr18:6348028\""] #<-- all inputs per bac file has to be formatted like this

#this should be the final output per bac file, each element in string is a whole string (not exceeding ) containing chrIDs
#of same chromosome for 1 -22 of chromosomes
#print(y)

#first make folder where csv Fies will go
    
folderName = makeAndGetFolderNameForBac(csv1)
parent_dir = "/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned"
path = os.path.join(parent_dir, folderName)
#os.mkdir(path) <--- file exists already -make exception

#new command every ~1000 chrIDs
    

count = 1
start = time.time()
for i in listOfStringUnixInputs:
    #print(i)
    selectEleToGetChrom = i.split("\\n")[1]
    chrom = selectEleToGetChrom.split(":")[0]
    filename = "%s_%s.csv"%(chrom,count)
    #print(filename)
    cmd = "curl -k -H \"Content-Type: application/json\" -X POST -d '{\"snps\": %s, \"pop\": \"EUR\", \"r2_threshold\": \"0.3\", \"maf_threshold\": \"0.01\"}' 'https://ldlink.nci.nih.gov/LDlinkRest/snpclip?token=8a70a1b4adf8' >> /Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned/LDLink_Output/%s/%s"%(i, folderName, filename)
    count += 1
    os.system(cmd)
end = time.time()
print("Time taken is  %s seconds"%(end - start))







####Serializing entireStringOfIDs for later use as input
import json

file = open("/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned/output_of_1bac.txt", "w")
json.dump(entireStringOfIds, file)
file.close()"""


