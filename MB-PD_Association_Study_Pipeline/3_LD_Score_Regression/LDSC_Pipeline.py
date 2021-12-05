#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 19:20:39 2021

@author: rodrigosandon
"""

# only choose bac csv files that end with .csv and that don't include "unknown" string
import pandas as pd
import csv
import os
import time
import os.path
from os import path

# the PD file for LD SR run is : /Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/nalls_onlyRsIDs.sumstats.gz

# Ex: /Volumes/Passport/119Bacs_addedP/addedPgenus..Clostridiuminnocuumgroup.id.14397.summary.txt.csv


def formatSumStatsForBac(csvPath):
    newFileName = "formatted." + csvPath.split("/")[4].replace(
        ".txt", "").replace("..", ".").replace("addedPgenus", "addedP.genus")
    newFilePath = "/Volumes/Passport/formatted119Bacs_addedP/%s" % (
        newFileName)

    if path.exists(newFilePath) == False:
        df = pd.read_csv(csvPath)
        newDf = df[['rsID', 'eff.allele', 'ref.allele',
                    'beta', 'SE', 'N', 'pDerived']].copy()
        newDf.insert(3, 'Zscore', newDf['beta'] / newDf['SE'], True)
        newDf = newDf.drop(['beta', 'SE'], axis=1)
        newDf = newDf.rename(columns={
                             'rsID': 'snpid', 'pDerived': 'P-value', 'eff.allele': 'A1', 'ref.allele': 'A2'})

        # Now formatted.addedP.genus.Clostridiuminnocuumgroup.id.14397.summary.csv

        newDf.to_csv(newFilePath, index=False)
    else:
        print("File %s already exists" % (newFilePath))
    return newFilePath


def CSVtoTXT(csv_file, txtFile):
    with open(txtFile, "w") as my_output_file:
        with open(csv_file, "r") as my_input_file:
            [my_output_file.write(" ".join(row)+'\n')
             for row in csv.reader(my_input_file)]
        my_output_file.close()

# ex csv file now: /Volumes/Passport/formatted119Bacs_addedP/formatted.addedP.genus.Clostridiuminnocuumgroup.id.14397.summary.csv


def mungeDataCall(csv_file):

    txtFilePath = csv_file.replace(".csv", ".txt")
    # /Volumes/Passport/formatted119Bacs_addedP/formatted.addedP.genus.Clostridiuminnocuumgroup.id.14397.summary.txt

    newFilePath = txtFilePath.replace(
        "formatted119Bacs_addedP", "munge119Bacs_output").replace(".txt", "")

    if path.exists(txtFilePath) == False:
        CSVtoTXT(csv_file, txtFilePath)
    else:
        print("File %s already munged" % (newFilePath))
    # os.chdir("/Users/rodrigosandon/ldsc")

    cmd = "./munge_sumstats.py \
        --sumstats %s \
        --out %s \
        --merge-alleles /Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/ldsc/w_hm3.snplist" % (txtFilePath, newFilePath)

    if path.exists(newFilePath + ".sumstats.gz") == False:  # only if .gz file don't exist
        os.system(cmd)

    return newFilePath

# ex gx bac name munged: /Volumes/Passport/munge119Bacs_output/formatted.addedP.genus.Clostridiuminnocuumgroup.id.14397.summary.sumstats.gz


def LDscore_regression(munged_bac_output):
    print("munged_bac_output: ", munged_bac_output)
    LDSR_outName = "/Volumes/T7Touch/NIHSummer2021/Code/LDSC_analysis2/results/%s_ldscResults" % (
        munged_bac_output.split("/")[4].split(".")[3])  # <--only bac name

    # os.chdir("/Users/rodrigosandon/ldsc")

    cmd = "./ldsc.py \
        --rg %s,/Volumes/T7Touch/NIHSummer2021/Code/LDSC_analysis2/munged_META5_all_with_rsid.sumstats.gz \
        --ref-ld-chr /Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/ldsc/eur_w_ld_chr/ \
        --w-ld-chr /Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/ldsc/eur_w_ld_chr/ \
        --out %s " % (munged_bac_output + ".sumstats.gz", LDSR_outName)

    os.system(cmd)

    return LDSR_outName

###MAIN###


masterDir = "/Volumes/Passport/119Bacs_addedP/"


for root, dirs, files in os.walk(masterDir):
    for name in files:
        start = time.time()
        bacPathToProcess = os.path.join(root, name)
        print("Processing", bacPathToProcess)

        newFilePath1 = formatSumStatsForBac(bacPathToProcess)

        newFilePath2 = mungeDataCall(newFilePath1)

        LDSR_outName = LDscore_regression(newFilePath2)

        print("Results can be found in", LDSR_outName)
        end = time.time()
        print("Time to perform LDSR:", (end-start)/60, "mins")

# BEFORE RUNNING FILE, MUNGE NEW PD GWAS w/ w_hm3 snplist
cmd = "./munge_sumstats.py \
        --sumstats /Volumes/T7Touch/NIHSummer2021/Code/LDSC_analysis2/reformatted2_META5_all_with_rsid.txt \
        --out /Volumes/T7Touch/NIHSummer2021/Code/LDSC_analysis2/munged_META5_all_with_rsid \
        --merge-alleles /Volumes/T7Touch/NIHSummer2021/Data/LDSR_analysis/ldsc/w_hm3.snplist"

# os.system(cmd)
# RUN THIS FILE ON TERMINAL:
# make sure the pd gwas is gziped
# make sure cd /Users/rodrigosandon/ldsc
# conda activate ldsc
# python --version should be python2
# python /Volumes/T7Touch/NIHSummer2021/Code/LDSC_Pipeline.py >> /Volumes/T7Touch/NIHSummer2021/Code/LDSC_analysis2/results/results.txt
