# -*- coding: utf-8 -*-
"""
2022.10.28 Heesoo SONG

Match spot positions
"""

import pandas as pd
import math
from collections import Counter



def eucledean (x1, x2, y1, y2):
    distance = math.sqrt(((x1 - x2)**2) + ((y1 - y2)**2))
    return distance

def Find_matches(df1, df2, ch1, ch2):
    match_dic = {}
    for i in range(0, df1.shape[0]):
        df1_x = df1.loc[i, 'x [nm]']
        df1_y = df1.loc[i, 'y [nm]']
        
        # Filter possible spots before calculating distance
        # Reduce number of calculations
        df_search = df2[(df2['x [nm]'] < (df1_x + limit)) & (df2['x [nm]'] > (df1_x - limit)) & 
                        (df2['y [nm]'] < (df1_y + limit)) & (df2['y [nm]'] > (df1_y - limit))]
        
        # Calculate distance between two spots
        for j in df_search.index:
            search_x = df_search.loc[j, 'x [nm]']
            search_y = df_search.loc[j, 'y [nm]']
            
            distance = eucledean(df1_x, search_x, df1_y, search_y)
            
            # Add to dictionary
            if distance <= limit:
                if i in list(match_dic.keys()):
                   item = match_dic.get(i)
                   item = item.append(j)
                else:
                    match_dic[i] = [j]
                    
                    
        
    # Check df1 spots that corresponded to multiple df2 spots
    # (one df1 -> multiple df2 spots)
    match_counts = {}
    for i in match_dic.keys():
        match_counts[i] = len(match_dic.get(i))
    
    multiple_matches = {k: v for k, v in match_counts.items() if v > 1}
    
    if len(multiple_matches) > 0:
        print("Multiple Matches (df1 -> multiple df2)")
        print(f"----- {multiple_matches}")
        
        # Compare distances between possible spots and find out closest spot
        # Then, update dictionary to only leave single match
        for i in multiple_matches.keys():
            
            # Find minimum distance
            compare_distance = []
            df2_matches = match_dic.get(i)
            for j in df2_matches:
                x1 = df1.loc[i, 'x [nm]']
                y1 = df1.loc[i, 'y [nm]']
                x2 = df2.loc[j, 'x [nm]']
                y2 = df2.loc[j, 'y [nm]']
                distance = eucledean(x1, x2, y1, y2)
                compare_distance.append(distance)
            min_index = compare_distance.index(min(compare_distance))
            
            # Update dictionary (remove multiple matches)
            match_dic[i] = [list(match_dic.get(i))[min_index]]
        print("----- Dictionary updated to avoid multiple matches \n")
    
    
    
    # Check if the spot matching is unique (multiple df1 -> one df2 spot)
    matched_spots_list = []
    for i in match_dic.keys():
        elements = match_dic.get(i)
        for j in elements:
            matched_spots_list.append(j)
            
    matched_spots_unique = set(matched_spots_list)
    
    if len(matched_spots_list) == len(matched_spots_unique):
        print("Unique matches (single df1 -> df2)")
    else:
        print("Multiple matches (multiple df1 -> df2):")
        
        # How many times the spot from df2 matched with the spot from df1
        spot_counts = Counter(matched_spots_list)
        non_unique_spots = {k: v for k, v in spot_counts.items() if v > 1}
        print(f"----- Number of occurrence {non_unique_spots}")
        
        # Which spots from df1 matched multiple times with the spot from df2
        for i in non_unique_spots.keys():
            multiple_matches = {k: v for k, v in match_dic.items() if i in v}
            print(f"---------- {multiple_matches}")
        
            # Compare distances between possible spots and find out closest spot
            # Then, update dictionary to only leave single match
            # Different from (df1 -> multiple df2), we will remove key:value pair
            min_distance = (limit + 1)
            min_key = 0
            for a, b in multiple_matches.items():
                # Find minimum distance
                x1 = df1.loc[a, 'x [nm]']
                y1 = df1.loc[a, 'y [nm]']
                x2 = df2.loc[b, 'x [nm]']
                y2 = df2.loc[b, 'y [nm]']
                distance = eucledean(x1, x2, y1, y2)
                if distance < min_distance:
                    min_distance = distance
                    min_key = a
    
            # Update dictionary (remove multiple matches)
            keys_to_remove = list(multiple_matches.keys())
            keys_to_remove.remove(min_key)
            for j in keys_to_remove:
                match_dic.pop(j)
        print("----- Dictionary updated to avoid multiple matches \n")
    
    return match_dic
        
        
# =============================================================================
# Integrate two dataframes
# =============================================================================

def add_extra_info(df_integrated, df_indices, df, ch):
    ## pos
    df_integrated[ch+'_posX'] = list(df.loc[df_indices, 'x [nm]'])
    df_integrated[ch+'_posY'] = list(df.loc[df_indices, 'y [nm]'])
    
    ## id
    df_integrated[ch+'_id'] = [x + 1 for x in df_indices]
    
    ## sigma
    df_integrated[ch+'_sigma'] = list(df.loc[df_indices, 'sigma [nm]'])
    
    ## intensity
    df_integrated[ch+'_int'] = list(df.loc[df_indices, 'intensity [photon]'])
    
    return df_integrated

def integrate_unmatched_spot_info(df_integrated, df_unmatched_indices, df, ch):
    temp_df = pd.DataFrame(columns=['x [nm]', 'y [nm]', 
                                          'C1_id', 'C2_id', 'C3_id', 'C4_id',
                                          'C1_posX', 'C1_posY', 'C2_posX', 'C2_posY', 
                                          'C3_posX', 'C3_posY', 'C4_posX', 'C4_posY',
                                          'C1_sigma', 'C2_sigma', 'C3_sigma', 'C4_sigma',
                                          'C1_int', 'C2_int', 'C3_int', 'C4_int'])
    
    temp_df['x [nm]'] = list(df.loc[df_unmatched_indices, 'x [nm]'])
    temp_df['y [nm]'] = list(df.loc[df_unmatched_indices, 'y [nm]'])
    
    temp_df = add_extra_info(temp_df, df_unmatched_indices, df, ch)
    
    df_integrated = pd.concat([df_integrated, temp_df])
    
    return df_integrated



base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"

df_C1 = pd.read_csv(base + "C1_Result.csv")
df_C2 = pd.read_csv(base + "C2_Result.csv")
df_C3 = pd.read_csv(base + "C3_Result.csv")
df_C4 = pd.read_csv(base + "C4_Result.csv")

limit = 150
#print(df1.columns)

## Create an empty DataFrame with column name
df_integrated = pd.DataFrame(columns=['x [nm]', 'y [nm]', 
                                      'C1_id', 'C2_id', 'C3_id', 'C4_id',
                                      'C1_posX', 'C1_posY', 'C2_posX', 'C2_posY', 
                                      'C3_posX', 'C3_posY', 'C4_posX', 'C4_posY',
                                      'C1_sigma', 'C2_sigma', 'C3_sigma', 'C4_sigma',
                                      'C1_int', 'C2_int', 'C3_int', 'C4_int'])
dataframes = [df_C2, df_C3, df_C4]
channels = ['C2', 'C3', 'C4']

df1 = df_C1
ch1 = 'C1'
for k in range(1,len(dataframes)):
    df2 = dataframes[k]
    ch2 = channels[k]
    
    # Find matches
    match_dic = Find_matches(df1, df2, ch1, ch2)


    # Integrate matching spots
    
    print(len(df_integrated))
    
    ## Average position
    df1_match_indices = list(match_dic.keys())
    df2_match_indices = [index for sublist in list(match_dic.values()) for index in sublist]
    
    match_x_df1 = list(df1.loc[df1_match_indices, 'x [nm]'])
    match_x_df2 = list(df2.loc[df2_match_indices, 'x [nm]'])
    
    match_y_df1 = list(df1.loc[df1_match_indices, 'y [nm]'])
    match_y_df2 = list(df2.loc[df2_match_indices, 'y [nm]'])
    
    average_x = []
    for i in range(0, len(match_x_df1)):
        avg = (match_x_df1[i] + match_x_df2[i])/2
        average_x.append(avg)
        
    average_y = []
    for i in range(0, len(match_y_df1)):
        avg = (match_y_df1[i] + match_y_df2[i])/2
        average_y.append(avg)
        
    df_integrated['x [nm]'] = average_x
    df_integrated['y [nm]'] = average_y
    
    
    ## Add other info
    df_integrated = add_extra_info(df_integrated, df1_match_indices, df1, ch1)
    df_integrated = add_extra_info(df_integrated, df2_match_indices, df2, ch2)
    
    
    # Integrate unmatching spots
    ## Average position = pos
    df1_unmatch_indices = [x for x in list(df1.index) if x not in list(match_dic.keys())]
    df2_unmatch_indices = [x for x in list(df2.index) if x not in df2_match_indices]
    
    df_integrated = integrate_unmatched_spot_info(df_integrated, df1_unmatch_indices, df1, ch1)
    df_integrated = integrate_unmatched_spot_info(df_integrated, df2_unmatch_indices, df2, ch2)

    break

