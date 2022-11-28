# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 14:46:49 2022

@author: user
"""

import pandas as pd
import math
import numpy

def eucledean (x1, x2, y1, y2):
    distance = math.sqrt(((x1 - x2)**2) + ((y1 - y2)**2))
    return distance

def set_to_merge(df, merge_limit, check_list, merge_set=[], merge_set_list = []):
    
    # Examine every spots for neighboring spots within merge_limit
    while len(check_list) > 0:
        df_x = df.loc[check_list[0], 'x [nm]']
        df_y = df.loc[check_list[0], 'y [nm]']
        
        # Filter possible spots before calculating distance
        # Reduce number of calculations
        df_search = df[(df['x [nm]'] < (df_x + merge_limit)) & (df['x [nm]'] > (df_x - merge_limit)) & 
                        (df['y [nm]'] < (df_y + merge_limit)) & (df['y [nm]'] > (df_y - merge_limit))]
        
        df_search = df_search.loc[df_search.index >= check_list[0]]
        
        # If there are any other spots within the range, evaluate the distance
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
                    
                    # Search if there is another neighboring spot of the neighboring spot (function within function).
                    # Since we only need to search for the neighboring spots of one spot (df_search.index[i]),
                    #   the length of check_list is 1.
                    merge_set, merge_set_list = set_to_merge(df, merge_limit, [df_search.index[i]], merge_set, merge_set_list)
        
        # Make sure this block is not operated in function within function, which we can distinguish by the length of check_list.
        # When examination of one spot is done, remove it from the check_list, so that we can end the while loop.
        if len(check_list) != 1:
            if len(merge_set) > 0:
                for spot in set(merge_set):
                    check_list.remove(spot)
                merge_set_list.append(set(merge_set))
                merge_set = []
                
            else:
                check_list.remove(check_list[0])
        else:
            check_list.remove(check_list[0])
    
    return merge_set, merge_set_list

def merge_multiple_spots_in_one_cell(df, merge_limit, merge_variable = "intensity/sigma", merge_option=min):
    check_list = list(range(df.shape[0]))
    merge_set, merge_set_list = set_to_merge(df, merge_limit, check_list)
    
    print(f"{len(merge_set_list)} spots with multiple marks")
    print(f"- Before merge = {df.shape}")
    
    # For all the merge_set lists, select one spot that can represent the others
    #   according to the options you have set.
    n_merge = 0
    for spots in merge_set_list:
        spots = list(spots)
        compare_list = []
        for x in spots:
            if merge_variable == "intensity [nm]" or merge_variable == "sigma [nm]":
                value = df.loc[x, merge_variable]
                compare_list.append(value)
                
            elif merge_variable == "intensity/sigma":
                intensity = df.loc[x, 'intensity [photon]']
                sigma = df.loc[x, 'sigma [nm]']
                
                compare_list.append(intensity/sigma)
                
        n_merge += (len(compare_list) - 1)
        
        # Replace x,y position by the average value of the spots to merge
        avg_x = numpy.mean(df.loc[spots, "x [nm]"])
        avg_y = numpy.mean(df.loc[spots, "y [nm]"])
        
        spot_to_remain = compare_list.index(merge_option(compare_list))
        
        df.loc[spots[spot_to_remain], "x [nm]"] = avg_x
        df.loc[spots[spot_to_remain], "y [nm]"] = avg_y
        
        # Drop the rows except the spot to maintain
        spots.pop(spot_to_remain)
        df.drop(spots, inplace=True)
    

    df.reset_index(drop=True, inplace=True)
    print(f"- After merge = {df.shape}    : {n_merge} spot info removed/merged \n")
    
    return df

if __name__ == "__main__":

    #base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"
    base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"
    #base = "C:/Users/user/Desktop/20221123_Pos10_Neg10/Positive/PB566_01/"
    
    df_C1 = pd.read_csv(base + "C1_Result.csv")
    df_C2 = pd.read_csv(base + "C2_Result.csv")
    df_C3 = pd.read_csv(base + "C3_Result.csv")
    df_C4 = pd.read_csv(base + "C4_Result.csv")
    
    # Three variables and two options to merge spots:
    #     "sigma [nm]", "intensity [nm]", "intensity/sigma"
    #     min, max
    # "intensity/sigma" and "min" as default
    df_C1 = merge_multiple_spots_in_one_cell(df_C4, 60, "intensity/sigma", min)