#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:56:49 2021

@author: rodrigosandon
"""

import pandas as pd
import gzip
import time
import os

gzPath = "/Volumes/Passport/MiBioGen_QmbQTL_summary_genus/genus..Clostridiuminnocuumgroup.id.14397.summary.txt.gz"

a_file = gzip.open(gzPath, "rb")
contents = a_file.read()

txtin = open(contents, "rt")
csvout = open("/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/MBsumStats.csv", "wt")

for line in txtin:
    csvout.write(' '.join(line.split()))
    csvout.write(' \n')

txtin.close()
csvout.close()


readFile = pd.read_csv(contents, 
                       delim_whitespace=True,  error_bad_lines=False)

readFile.to_csv("/Volumes/T7Touch/NIHSummer2021/Data/MicrobiomeVariants/MBsumStats.csv", index=None)