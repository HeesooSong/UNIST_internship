# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd

base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Positive/3/Results/"

C1_df = pd.read_csv(base + "C1_raw_64nm_Result.csv")
C1_df

C3_df = pd.read_csv(base + "C3_raw_64nm_Result.csv")