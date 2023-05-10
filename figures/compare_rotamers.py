#!/usr/bin/env python3

"""

"""
from qfit import Structure
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import warnings
warnings.filterwarnings("ignore")

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
        return diff < 35


    TP_count, loss_count, gain_count, FP_count, match_count, no_match, shared_different = 0, 0, 0, 0, 0, 0, 0
   #print(df1)
   # print(df2)

    # Add a new column 'matched' to df1 and df2 and set all values to False
    df1['matched'] = False
    df2['matched'] = False
    for index1, row1_single in df1.iterrows():
            if row1_single['matched']:  # Skip rows that have already been matched
                #print('continue')
                continue
            chi1_1, chi2_1 = row1_single['chi1'], row1_single['chi2']
        
            for index2, row2_single in df2.iterrows():
                #print(index2)
                if row2_single['matched']:  # Skip rows that have already been matched
                    #print('continue')
                    continue
                chi1_2, chi2_2 = row2_single['chi1'], row2_single['chi2']

                # Check if the angles of row1 and row2 match
                if pd.isnull(chi2_1) or pd.isnull(chi2_2):
                    if angle_difference_greater_than_15(chi1_1, chi1_2):
                        match_count = 1
                        #print('match')
                        df1.loc[index1, 'matched'] = True
                        df2.loc[index2, 'matched'] = True
                        # If there's a match and more than one row1_single, move on to the next row1_single
                    else:
                        no_match = 1
                       #print('no match')
                else:                  
                    if (angle_difference_greater_than_15(chi1_1, chi1_2) and angle_difference_greater_than_15(chi2_1, chi2_2)):
                        match_count = 1
                        df1.loc[index1, 'matched'] = True
                        df2.loc[index2, 'matched'] = True
                        #print('match')
                        # If there's a match and more than one row1_single, move on to the next row1_single
                    else:
                        no_match = 1
                        #print('no match')
        
        # Save the evaluation results
    if no_match == 0:
        TP_count = 1
    elif match_count == 0:
        FP_count = 1
        print(df1)
        print(df2)
    elif no_match > 0 and match_count > 0:
        if len(df1) > len(df2):
            loss_count = 1
        elif len(df2) > len(df1):
            gain_count = 1
        else:
            TP_count = 1
    else:
        print('ambigous')

    # Return the evaluation results
    return {'TP': TP_count, 'loss': loss_count, 'gain': gain_count, 'FP': FP_count, 'share_diff': shared_different}

def plot_results(results, alt):
    def calc_sum(results_df):
            # Calculate the sum of TP, FP, Gain, and Loss columns
            total_count = results_df["TP"].sum() + results_df["FP"].sum() + results_df["gain"].sum() + results_df["loss"].sum() + results_df["share_diff"].sum()
            # Calculate the percentages of TP, FP, Gain, and Loss
            results_df["TP_pct"] = results_df["TP"] / total_count * 100
            results_df["FP_pct"] = results_df["FP"] / total_count * 100
            results_df["gain_pct"] = results_df["gain"] / total_count * 100
            results_df["loss_pct"] = results_df["loss"] / total_count * 100
            results_df["share_diff"] = results_df["share_diff"] / total_count * 100
            # Print the percentages
            print("TP %:", results_df["TP_pct"].sum())
            print("FP %:", results_df["FP_pct"].sum())
            print("Gain %:", results_df["gain_pct"].sum())
            print("Loss %:", results_df["loss_pct"].sum())
            print("Share Some %:", results_df["share_diff"].sum())

            all_array = [results_df["TP_pct"].sum(), results_df["FP_pct"].sum(), results_df["gain_pct"].sum(),results_df["loss_pct"].sum(), results_df["share_diff"].sum()]
            return all_array



    def save_legend():
        # Create a dummy figure with the same colors and labels as the main plot
    # Create a dummy figure with the same colors and labels as the main plot
        fig_legend, ax_dummy = plt.subplots()
        ax_dummy.set_visible(False)

        # Add bars and labels to the dummy figure
        labels = ['True Positives', 'False Positives', 'Gain', 'Loss', 'Share Some']
        colors = ['#191970', '#097969', '#FF5733', '#C70039', '#FFC300']

        # Create custom patches for the legend
        legend_elements = [Patch(facecolor=color, label=label) for color, label in zip(colors, labels)]

        # Create the legend using custom patches
        legend = plt.legend(handles=legend_elements, bbox_to_anchor=(0.5, 0.5), loc='center', frameon=False)

        # Save the legend to a separate file
        fig_legend.canvas.draw()
        fig_legend.savefig('RotamerStatus_legend.png', bbox_inches='tight', dpi=300)
        plt.close(fig_legend)



    print('All:')
    all_ = calc_sum(results)
    print('Alt:')
    alt_ = calc_sum(alt)  

    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)

    # Create a stacked bar plot for the first set of values
        # Adjust the width of the bars and the space between them
    bar_width = 0.03
    space_between = 0.1

    # Create a stacked bar plot for the first set of values
    all_bar = ax1.bar('All', all_[0], color='#191970', width=bar_width)
    all_bar = ax1.bar('All', all_[1], bottom=all_[0], color='#097969', width=bar_width)
    all_bar = ax1.bar('All', all_[2], bottom=all_[0] + all_[1], color='#FF5733', width=bar_width)
    all_bar = ax1.bar('All', all_[3], bottom=all_[0] + all_[1] + all_[2], color='#C70039', width=bar_width)
    all_bar = ax1.bar('All', all_[4], bottom=all_[0] + all_[1] + all_[2] + all_[3], color='#FFC300', width=bar_width)

    # Create a stacked bar plot for the second set of values
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[0], color='#191970', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[1], bottom=alt_[0], color='#097969', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[2], bottom=alt_[0] + alt_[1], color='#FF5733', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[3], bottom=alt_[0] + alt_[1] + alt_[2], color='#C70039', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[4], bottom=alt_[0] + alt_[1] + alt_[2] + alt_[3], color='#FFC300', width=bar_width)

    # Call the save_legend function between the plot creation and display


    # Adjust the space between the two bar plots
    plt.subplots_adjust(wspace=space_between)

    # Add a common y-axis label
    fig.text(0.04, 0.5, 'Percentage', va='center', rotation='vertical', fontsize=12)

    # Save the bar plot to a file
    plt.savefig(f'RotamerStatus_stacked_bars.png')
    save_legend()



def main():
    qfit = read_files('/Users/stephaniewanko/Downloads/temp/qfit_test_set', 'nconfs', 58,62)
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
    print(results_df.head())

    # Save the DataFrame to a CSV file
    results_df.to_csv("comparison_results.csv", index=False)
    results_df_alt.to_csv("comparison_results_alt.csv", index=False)


    # Call the plot_results function
    plot_results(results_df, results_df_alt)


        

if __name__ == "__main__":
    main()
