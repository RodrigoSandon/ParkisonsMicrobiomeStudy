#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 19:20:39 2021

@author: rodrigosandon
"""

import pandas as pd

def formatSumStatsForLDSC_PD (csvPath):
    df = pd.read_csv(csvPath)
    newDf = df[['SNP', 'A1', 'A2', 'b','se', 'N_cases', 'p']].copy()
    newDf.insert(3, 'Zscore', newDf['b'] / newDf['se'], True)
    newDf = newDf.drop(['b','se'], axis = 1)
    newDf = newDf.rename(columns={'N_cases' : 'N', 'SNP' : 'snpid', 'p' : 'P-value'}) # Now have: snpid A1 A2 Zscore N P-value
    
def formatSumStatsForLDSC_MB (csvPath):
    df = pd.read_csv(csvPath)
    
    
chrID_rsIDs = pd.read_csv("/Volumes/T7Touch/NIHSummer2021/Data/EuropeanParkinsons/chrPosRs.tab", sep='\t', lineterminator='\r', header = None)
chrID_rsIDs.head()