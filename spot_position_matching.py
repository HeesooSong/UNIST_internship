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

def Find_matches(df1, df2):
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
        print("----- Multiple Matches (df1 -> multiple df2)")
        print(f"---------- {multiple_matches}")
        
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
        print("---------- Dictionary updated to avoid multiple matches \n")
    else:
        print("----- Unique Matches (df1 -> single df2")
      
    
    # Check if the spot matching is unique (multiple df1 -> one df2 spot)
    matched_spots_list = []
    for i in match_dic.keys():
        elements = match_dic.get(i)
        for j in elements:
            matched_spots_list.append(j)
            
    matched_spots_unique = set(matched_spots_list)
    
    if len(matched_spots_list) == len(matched_spots_unique):
        print("----- Unique matches (single df1 -> df2)")
    else:
        print("----- Multiple matches (multiple df1 -> df2):")
        
        # How many times the spot from df2 matched with the spot from df1
        spot_counts = Counter(matched_spots_list)
        non_unique_spots = {k: v for k, v in spot_counts.items() if v > 1}
        print(f"---------- Number of occurrence {non_unique_spots}")
        
        # Which spots from df1 matched multiple times with the spot from df2
        for i in non_unique_spots.keys():
            multiple_matches = {k: v for k, v in match_dic.items() if i in v}
            print(f"--------------- {multiple_matches}")
        
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
        print("---------- Dictionary updated to avoid multiple matches \n")
    
    return match_dic
        

def add_extra_info(df_integrated, integrated_indices, df_indices, df, ch):
    ## pos
    df_integrated.loc[integrated_indices, ch+'_posX'] = list(df.loc[df_indices, 'x [nm]'])
    df_integrated.loc[integrated_indices, ch+'_posY'] = list(df.loc[df_indices, 'y [nm]'])
    
    ## id
    df_integrated.loc[integrated_indices, ch+'_id'] = [x + 1 for x in df_indices]
    
    ## sigma
    df_integrated.loc[integrated_indices, ch+'_sigma'] = list(df.loc[df_indices, 'sigma [nm]'])
    
    ## intensity
    df_integrated.loc[integrated_indices, ch+'_int'] = list(df.loc[df_indices, 'intensity [photon]'])
    
    return df_integrated


def integrate_unmatched_spot_info(df_integrated, df_unmatched_indices, df, ch):
    temp_df = pd.DataFrame(columns=['x [nm]', 'y [nm]', 'n_match',
                                'C1_id', 'C2_id', 'C3_id', 'C4_id',
                                'C1_posX', 'C1_posY', 'C2_posX', 'C2_posY', 
                                'C3_posX', 'C3_posY', 'C4_posX', 'C4_posY',
                                'C1_sigma', 'C2_sigma', 'C3_sigma', 'C4_sigma',
                                'C1_int', 'C2_int', 'C3_int', 'C4_int'])
    
    temp_df['x [nm]'] = list(df.loc[df_unmatched_indices, 'x [nm]'])
    temp_df['y [nm]'] = list(df.loc[df_unmatched_indices, 'y [nm]'])
    
    temp_indices = range(len(df_unmatched_indices))
    
    temp_df = add_extra_info(temp_df, temp_indices, df_unmatched_indices, df, ch)
    
    df_integrated = pd.concat([df_integrated, temp_df])
    
    return df_integrated

def average_position(df_integrated, indices):
    average_tempX = df_integrated[['C1_posX', 'C2_posX', 'C3_posX', 'C4_posX']]
    average_tempY = df_integrated[['C1_posY', 'C2_posY', 'C3_posY', 'C4_posY']]
    
    df_integrated.loc[indices, 'x [nm]'] = average_tempX.mean(axis=1)
    df_integrated.loc[indices, 'y [nm]'] = average_tempY.mean(axis=1)
    
    return df_integrated


def initialize_df_integrated(df, ch):
    # Create an empty DataFrame with column name
    df_integrated = pd.DataFrame(columns=['x [nm]', 'y [nm]', 'n_match',
                                          'C1_id', 'C2_id', 'C3_id', 'C4_id',
                                          'C1_posX', 'C1_posY', 'C2_posX', 'C2_posY', 
                                          'C3_posX', 'C3_posY', 'C4_posX', 'C4_posY',
                                          'C1_sigma', 'C2_sigma', 'C3_sigma', 'C4_sigma',
                                          'C1_int', 'C2_int', 'C3_int', 'C4_int'])
    # Add initial info
    df_integrated['x [nm]'] = df['x [nm]']
    df_integrated['y [nm]'] = df['y [nm]']
    
    df_indices = range(len(df_integrated))
    
    df_integrated = add_extra_info(df_integrated, df_indices, df_indices, df, ch)
    
    return df_integrated


def main(df_C1, df_C2, df_C3, df_C4, limit):
    
    dataframes = [df_C2, df_C3, df_C4]
    channels = ['C2', 'C3', 'C4']
    
    ch1 = 'C1'
    
    df_integrated = initialize_df_integrated(df_C1, ch1)
    df1 = df_integrated
    
    for k in range(0,len(dataframes)):
        df2 = dataframes[k]
        ch2 = channels[k]
        
        print(f'Integrating {ch1} and {ch2}')
        print(f'----- Size of {ch1} / {ch2}: {len(df1)} / {len(df2)}')
        
        # Find matches
        match_dic = Find_matches(df1, df2)
        
        # Integrate matching spots    
        ## Average position
        df1_match_indices = list(match_dic.keys())
        df2_match_indices = [index for sublist in list(match_dic.values()) for index in sublist]

            
        print(f'----- Matched positions: {len(df1_match_indices)}')
        print(f'----- Size bf integration: {len(df_integrated)}')

        # Add other info
        df_integrated = add_extra_info(df_integrated, df1_match_indices, df2_match_indices, df2, ch2)
        
        # Average position
        df_integrated = average_position(df_integrated, df1_match_indices)

        # Integrate unmatching spots
        ## Average position = pos
        df2_unmatch_indices = [x for x in list(df2.index) if x not in df2_match_indices]
        
        df_integrated = integrate_unmatched_spot_info(df_integrated, df2_unmatch_indices, df2, ch2)
        
        print(f'----- Size after integration: {len(df_integrated)} \n')
    
        # Reset index
        df_integrated = df_integrated.reset_index(drop=True)
        
        # Set df1
        df1 = df_integrated
        ch1 = 'df_integrated'
    
    # Sort rows with number of matching
    df_integrated['n_match'] = df_integrated[['C1_id', 'C2_id', 'C3_id', 'C4_id']].count(axis=1)
    df_integrated = df_integrated.sort_values(by=['n_match'], ascending=False)
    
    # Reset index
    df_integrated = df_integrated.reset_index(drop=True)
    
    return df_integrated

def fill_missing_values(df_integrated, df_C1, df_C2, df_C3, df_C4):
    # Fill intensity values with average background signal intensity
    bkg_C1 = df_C1['bkgstd [photon]'].mean()
    bkg_C2 = df_C2['bkgstd [photon]'].mean()
    bkg_C3 = df_C3['bkgstd [photon]'].mean()
    bkg_C4 = df_C4['bkgstd [photon]'].mean()
    
    df_integrated['C1_int'] = df_integrated['C1_int'].fillna(bkg_C1)
    df_integrated['C2_int'] = df_integrated['C2_int'].fillna(bkg_C2)
    df_integrated['C3_int'] = df_integrated['C3_int'].fillna(bkg_C3)
    df_integrated['C4_int'] = df_integrated['C4_int'].fillna(bkg_C4)
    
    # Fill sigma with 0
    df_integrated[['C1_sigma', 'C2_sigma', 'C3_sigma', 'C4_sigma']] = df_integrated[['C1_sigma', 'C2_sigma', 'C3_sigma', 'C4_sigma']].fillna(0)
    
    
    
    return df_integrated
        
if __name__ == "__main__":

    base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"
    
    df_C1 = pd.read_csv(base + "C1_Result.csv")
    df_C2 = pd.read_csv(base + "C2_Result.csv")
    df_C3 = pd.read_csv(base + "C3_Result.csv")
    df_C4 = pd.read_csv(base + "C4_Result.csv")
    
    limit = 150
    
    df_integrated = main(df_C1, df_C2, df_C3, df_C4, limit)
    df_integrated.to_csv(base + "Spot_matching_result.csv")
    
    df_integrated = fill_missing_values(df_integrated, df_C1, df_C2, df_C3, df_C4)
    df_integrated.to_csv(base + "Spot_matching_result_imputated.csv")