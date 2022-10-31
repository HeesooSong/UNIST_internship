# -*- coding: utf-8 -*-
"""
2022.10.28 Heesoo SONG

Match spot positions
"""

import pandas as pd
import math
from collections import Counter


base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Positive/3/Results/"

df1 = pd.read_csv(base + "C1_raw_64nm_Result.csv")
df2 = pd.read_csv(base + "C3_raw_64nm_Result.csv")

limit = 150
#print(df1.columns)

def eucledean (x1, x2, y1, y2):
    distance = math.sqrt(((x1 - x2)**2) + ((y1 - y2)**2))
    return distance

# =============================================================================
# Find Matches between two dataframes
# =============================================================================
match_dic = {}
for i in range(1, df1.shape[0]):
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
    print(f"Multiple matches (multiple df1 -> df2):")
    
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
        
        
# =============================================================================
# Integrate two dataframes
# =============================================================================

average_position