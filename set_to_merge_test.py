# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 14:46:49 2022

@author: user
"""

import pandas as pd
import math

def eucledean (x1, x2, y1, y2):
    distance = math.sqrt(((x1 - x2)**2) + ((y1 - y2)**2))
    return distance

def set_to_merge(df, merge_limit, check_list, merge_set=[], merge_set_list = [], flag = False):
    while len(check_list) > 0:
        df_x = df.loc[check_list[0], 'x [nm]']
        df_y = df.loc[check_list[0], 'y [nm]']
        
        # Filter possible spots before calculating distance
        # Reduce number of calculations
        df_search = df[(df['x [nm]'] < (df_x + merge_limit)) & (df['x [nm]'] > (df_x - merge_limit)) & 
                        (df['y [nm]'] < (df_y + merge_limit)) & (df['y [nm]'] > (df_y - merge_limit))]
        
        df_search = df_search.loc[df_search.index >= check_list[0]]
        
        if len(df_search) > 1:
            
            spot1_x = df_search.loc[df_search.index[0], 'x [nm]']
            spot1_y = df_search.loc[df_search.index[0], 'y [nm]']
            
            for i in range(1, len(df_search)):
                spot2_x = df_search.loc[df_search.index[i], 'x [nm]']
                spot2_y = df_search.loc[df_search.index[i], 'y [nm]']
                
                distance = eucledean(spot1_x, spot2_x, spot1_y, spot2_y)
                
                if distance <= merge_limit:
                    if merge_set == []:
                        merge_set = [df_search.index[0], df_search.index[i]]
                    else:
                        merge_set.append(df_search.index[i])
                    merge_set, merge_set_list = set_to_merge(df, merge_limit, list(range(df_search.index[i], df.shape[0])), merge_set)
                
                check_list.remove(df_search.index[i])

            break
        
        check_list.remove(check_list[0])
    
    return merge_set, merge_set_list

def merge_multiple_spots_in_one_cell(df, merge_limit):
    check_list = list(range(df.shape[0]))
    merge_set = set_to_merge(df, merge_limit, check_list)
        
    
    return df

if __name__ == "__main__":

    #base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"
    base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"
    #base = "C:/Users/user/Desktop/20221123_Pos10_Neg10/Positive/PB566_01/"
    
    df_C1 = pd.read_csv(base + "C1_Result.csv")
    df_C2 = pd.read_csv(base + "C2_Result.csv")
    df_C3 = pd.read_csv(base + "C3_Result.csv")
    df_C4 = pd.read_csv(base + "C4_Result.csv")
    
    merge_multiple_spots_in_one_cell(df_C1, 50)