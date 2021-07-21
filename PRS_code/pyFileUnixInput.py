#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 23:50:15 2021

@author: rodrigosandon
"""

import json

json_file = open("/Volumes/T7Touch/NIHSummer2021/Data/MiBioGen_pruned/output_of_1bac.txt")

variable = json.load(json_file)
json_file.close()

print(variable)