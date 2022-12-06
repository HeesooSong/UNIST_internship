# -*- coding: utf-8 -*-
"""
2022.10.28 Heesoo SONG

Match spot positions
"""

import pandas as pd
import math
from collections import Counter
import numpy
import sys
import datetime



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

def merge_multiple_spots_in_one_cell(df, merge_limit, merge_variable = "intensity/sigma", merge_option=max):
    check_list = list(range(df.shape[0]))
    merge_set, merge_set_list = set_to_merge(df, merge_limit, check_list, merge_set=[], merge_set_list = [])
    
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




def Find_matches(df1, df2, limit):
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
        #print(f"---------- {multiple_matches}")
        
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
        print("---------- Dictionary updated to avoid multiple matches")
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
            #print(f"--------------- {multiple_matches}")
        
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
        print("---------- Dictionary updated to avoid multiple matches")
    
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


def integrate_unmatched_spot_info(df_integrated, df_unmatched_indices, df, ch, channels):
    # Column names according to the length of channels
    column_id = [x + "_id" for x in channels]
    column_posX = [x + "_posX" for x in channels]
    column_posY = [x + "_posY" for x in channels]
    column_sigma = [x + "_sigma" for x in channels]
    column_int = [x + "_int" for x in channels]
    columns = ['x [nm]', 'y [nm]', 'n_match'] + column_id + column_posX + column_posY + column_sigma + column_int

    temp_df = pd.DataFrame(columns=columns)
    
    temp_df['x [nm]'] = list(df.loc[df_unmatched_indices, 'x [nm]'])
    temp_df['y [nm]'] = list(df.loc[df_unmatched_indices, 'y [nm]'])
    
    temp_indices = range(len(df_unmatched_indices))
    
    temp_df = add_extra_info(temp_df, temp_indices, df_unmatched_indices, df, ch)
    
    df_integrated = pd.concat([df_integrated, temp_df])
    
    return df_integrated

def average_position(df_integrated, indices):
    n_channels = int((df_integrated.shape[1] - 3)/6)
    posX = []
    posY = []
    for i in range(1, n_channels+1):
        posX.append("C"+str(i)+"_posX")
        posY.append("C"+str(i)+"_posY")

    average_tempX = df_integrated[posX]
    average_tempY = df_integrated[posY]
    
    df_integrated.loc[indices, 'x [nm]'] = average_tempX.mean(axis=1)
    df_integrated.loc[indices, 'y [nm]'] = average_tempY.mean(axis=1)
    
    return df_integrated


def initialize_df_integrated(df, channels):
    # Column names according to the length of channels
    column_id = [x + "_id" for x in channels]
    column_posX = [x + "_posX" for x in channels]
    column_posY = [x + "_posY" for x in channels]
    column_sigma = [x + "_sigma" for x in channels]
    column_int = [x + "_int" for x in channels]
    columns = ['x [nm]', 'y [nm]', 'n_match'] + column_id + column_posX + column_posY + column_sigma + column_int
    
    # Create an empty DataFrame with column name
    df_integrated = pd.DataFrame(columns=columns)
    # Add initial info
    df_integrated['x [nm]'] = df['x [nm]']
    df_integrated['y [nm]'] = df['y [nm]']
    
    df_indices = range(len(df_integrated))
    
    df_integrated = add_extra_info(df_integrated, df_indices, df_indices, df, channels[0])
    
    return df_integrated


def main(df_list, channels, limit):

    df_integrated = initialize_df_integrated(df_list[0], channels)
    df1 = df_integrated
    ch1 = channels[0]
    
    for k in range(1,len(df_list)):
        df2 = df_list[k]
        ch2 = channels[k]
        
        print(f'Integrating {ch1} and {ch2}')
        print(f'----- Size of {ch1} / {ch2}: {len(df1)} / {len(df2)}')
        
        # Find matches
        match_dic = Find_matches(df1, df2, limit)
        
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
        
        df_integrated = integrate_unmatched_spot_info(df_integrated, df2_unmatch_indices, df2, ch2, channels)
        
        print(f'----- Size after integration: {len(df_integrated)} \n')
    
        # Reset index
        df_integrated = df_integrated.reset_index(drop=True)
        
        # Set df1
        df1 = df_integrated
        ch1 = 'df_integrated'

    channel_id = []
    for i in range(1, len(channels)+1):
        df_integrated['C'+str(i)+'_int_sig'] = df_integrated['C'+str(i)+'_int'] / df_integrated['C'+str(i)+'_sigma'] / df_integrated['C'+str(i)+'_sigma']
        channel_id.append('C'+str(i)+'_id')

    # Sort rows with number of matching
    df_integrated['n_match'] = df_integrated[channel_id].count(axis=1)
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

    base = "C:/Users/pc/Desktop/UNIST_internship/Sample_Image/Negative/2/"
    #base = "C:/Users/pc/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"

    df_C1 = pd.read_csv(base + "C1_Result.csv")
    df_C2 = pd.read_csv(base + "C2_Result.csv")
    df_C3 = pd.read_csv(base + "C3_Result.csv")
    df_C4 = pd.read_csv(base + "C4_Result.csv")

    # Text file where outputs exported
    #current_time = datetime.datetime.now()
    #date_time = current_time.strftime("%Y%m%d_%H%M%S")
    #sys.stdout = open(base + "log_" + date_time +".txt", "w")
    
    '''
    Merge the detection results if they are closer than the merge limit.
    In case of multiple detection, choose meta data of one spot that has min/max of sigma/intensity/(intensity/sigma)
    Three variables and two options to merge spots:
       - "sigma [nm]", "intensity [nm]", "intensity/sigma"(default)
       - min, max(default)
    '''
    merge_limit = 60
    df_C1 = merge_multiple_spots_in_one_cell(df_C1, merge_limit, "intensity/sigma", max)
    df_C2 = merge_multiple_spots_in_one_cell(df_C2, merge_limit, "intensity/sigma", max)
    df_C3 = merge_multiple_spots_in_one_cell(df_C3, merge_limit, "intensity/sigma", max)
    df_C4 = merge_multiple_spots_in_one_cell(df_C4, merge_limit, "intensity/sigma", max)

    df_C1.to_csv(base + "C1_Result_corrected.csv")
    df_C2.to_csv(base + "C2_Result_corrected.csv")
    df_C3.to_csv(base + "C3_Result_corrected.csv")
    df_C4.to_csv(base + "C4_Result_corrected.csv")

    # Perform spot matching based on the corrected/uncorrected spots
    df_list = [df_C1, df_C2, df_C3, df_C4] # You can define the order of integration here
    channels = ['C1', 'C2', 'C3', 'C4'] # Number of channels must match the length of df_list
    
    limit = 150
    
    df_integrated = main(df_list, channels, limit)
    df_integrated.to_csv(base + "Spot_matching_result.csv")

    df_integrated = fill_missing_values(df_integrated, df_C1, df_C2, df_C3, df_C4)
    df_integrated.to_csv(base + "Spot_matching_result_imputated.csv")

    #sys.stdout.close()