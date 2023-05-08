#!/usr/bin/env python3

"""

"""
from qfit import Structure
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt

def read_files(base_dir, comparison, file_start, file_end):
    all_files = glob.glob(base_dir + '/' + comparison + "/*_rotamer_output.txt")
    li = []
    for filename in all_files:
        try:
            df = pd.read_csv(filename, index_col=None, header=0, sep=':')
        except pd.errors.EmptyDataError:
            continue  
        df['PDB'] = filename[file_start:file_end]
        li.append(df)
    rotamer = pd.concat(li, axis=0, ignore_index=True)
    rotamer['category'] = comparison
    rotamer = rotamer[rotamer['residue']!= 'SUMMARY'].reset_index()
    split = rotamer['residue'].str.split(" ")
    for i in range(0, (len(rotamer.index)-1)):
        rotamer.loc[i,'chain'] = split[i][1]
        STUPID = str(rotamer.loc[i,'residue'])[3:6]
        tmp = []
        try:
            tmp = (int(''.join(i for i in STUPID if i.isdigit())))
        except:
            newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in STUPID)
            tmp = [float(i) for i in newstr.split()]
        rotamer.loc[i, 'resi'] = tmp
    return rotamer

def compare_rotamers(df1, df2):
    # Function to calculate the percentage difference between two angles
    def angle_difference(angle1, angle2):
        diff = abs(angle1 - angle2)
        return min(diff, 360 - diff)


# Function to check if angle1 and angle2 are different by more than 15 in either direction
    def angle_difference_greater_than_15(angle1, angle2):
        # If angle1 or angle2 is less than 180, subtract it from 360
        if angle1 < 180:
            angle1 = 360 - angle1
        if angle2 < 180:
            angle2 = 360 - angle2

        diff = abs(angle1 - angle2) 
        return diff < 30


    TP_count, loss_count, gain_count = 0, 0, 0

    # Iterate through each row in row1
    TP_count = 0
    FP_count = 0
    gain_count = 0
    loss_count = 0
    for _, row1_single in df1.iterrows():
        chi1_1, chi2_1 = row1_single['chi1'], row1_single['chi2']
        match_count, additional_row1, additional_row2 = 0, 0, 0
        
        # Iterate through each row in row2
        for _, row2_single in df2.iterrows():
            chi1_2, chi2_2 = row2_single['chi1'], row2_single['chi2']
            
            # Check if the angles of row1 and row2 match
            if pd.isnull(chi2_1) or pd.isnull(chi2_2):
                if angle_difference_greater_than_15(chi1_1, chi1_2):
                    match_count += 1
                    # If there's a match and more than one row1_single, move on to the next row1_single
                    if len(df1) > 1:
                        break
            else:                  
                if (angle_difference_greater_than_15(chi1_1, chi1_2) and angle_difference_greater_than_15(chi2_1, chi2_2)):
                    match_count += 1
                    if len(df1) > 1:
                        break
                else:
                    additional_row1 += 1
        
        # Save the evaluation results
    if match_count > 0:
        if additional_row1 == 0 and additional_row2 == 0:
            TP_count += 1
        elif additional_row1 > 0:
            gain_count += 1
        elif additional_row2 > 0:
            loss_count += 1
    else:
        FP_count += 1

    # Return the evaluation results
    return {'TP': TP_count, 'loss': loss_count, 'gain': gain_count, 'FP': FP_count}

def plot_results(results, alt):
    def calc_sum(results_df):
            # Calculate the sum of TP, FP, Gain, and Loss columns
            total_count = results_df["TP"].sum() + results_df["FP"].sum() + results_df["gain"].sum() + results_df["loss"].sum()

            # Calculate the percentages of TP, FP, Gain, and Loss
            results_df["TP_pct"] = results_df["TP"] / total_count * 100
            results_df["FP_pct"] = results_df["FP"] / total_count * 100
            results_df["gain_pct"] = results_df["gain"] / total_count * 100
            results_df["loss_pct"] = results_df["loss"] / total_count * 100
                    # Print the percentages
            print("TP %:", results_df["TP_pct"].sum())
            print("FP %:", results_df["FP_pct"].sum())
            print("Gain %:", results_df["gain_pct"].sum())
            print("Loss %:", results_df["loss_pct"].sum())

            all_array = [results_df["TP_pct"], results_df["FP_pct"], results_df["gain_pct"]]
            return all_array

    print('All:')
    all_ = calc_sum(results)
    print('Alt:')
    alt_ = calc_sum(alt)  

    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)

    # Create a stacked bar plot for the first set of values
    ax1.bar('All', all_[0], color='#191970')
    ax1.bar('All', all_[1], bottom=all_[0], color='#097969')
    ax1.bar('All', all_[2], bottom=all_[0] + all_[1], color='#FF5733')

    # Create a stacked bar plot for the second set of values
    ax2.bar('Deposited AltLocs', alt_[0], color='#191970', width=0.5)
    ax2.bar('Deposited AltLocs', alt_[1], bottom=alt_[0], color='#097969', width=0.5)
    ax2.bar('Deposited AltLocs', alt_[2], bottom=alt_[0] + altlocs_array[1], color='#FF5733', width=0.5)

    # Add a common y-axis label
    fig.text(0.04, 0.5, 'Percentage', va='center', rotation='vertical', fontsize=12)

    plt.savefig(f'RotamerStatus_stacked_half_width.png')




def main():
    qfit = read_files('/Users/stephaniewanko/Downloads/temp/qfit_test_set', 'final', 57,61)
    single = read_files('/Users/stephaniewanko/Downloads/temp/qfit_test_set', 'single_conf', 63, 67)
    unique_pdbs = set(single['PDB'].unique()).intersection(qfit['PDB'].unique())

    # Define an empty list to store comparison results
    results = []
    results_alt = []
    for pdb in unique_pdbs:
        # Filter the dataframes to only include rows with the current PDB
        df1_filtered = single[single['PDB'] == pdb]
        df2_filtered = qfit[qfit['PDB'] == pdb]

        # Iterate through each unique residue in the filtered dataframes
        unique_residues = set(df1_filtered['resi'].unique()).intersection(df2_filtered['resi'].unique())
        for resi in unique_residues:
            # Filter the dataframes to only include rows with the current residue
            df1_residue = df1_filtered[df1_filtered['resi'] == resi]
            df2_residue = df2_filtered[df2_filtered['resi'] == resi]

            comparison_result = compare_rotamers(df1_residue, df2_residue)
            comparison_result['residue'] = resi
            comparison_result['PDB'] = pdb
            # Check if df1_residue has more than 1 line
            if len(df1_residue) > 1:
            # Compare the rotamers of the current residue in both dataframes
                results.append(comparison_result)
                results_alt.append(comparison_result)
            else:
                results.append(comparison_result)

    # Convert results to a dataframe and print it
    results_df = pd.DataFrame(results)
    results_df_alt = pd.DataFrame(results_alt)

    # Save the DataFrame to a CSV file
    results_df.to_csv("comparison_results.csv", index=False)



    # Call the plot_results function
    plot_results(results_df, results_df_alt)


        

if __name__ == "__main__":
    main()
