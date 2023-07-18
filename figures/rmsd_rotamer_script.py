#!/usr/bin/env python3
"""

"""
#from qfit import Structure
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import warnings
warnings.filterwarnings("ignore")
plt.rcParams['figure.dpi'] = 600

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

    TP_count, loss_count, gain_count, FP_count, match_count, no_match, share = 0, 0, 0, 0, 0, 0, 0


    # Add a new column 'matched' to df1 and df2 and set all values to False
    df1['matched'] = False
    df2['matched'] = False
    continue_outer_loop = False
    for index1, row1_single in df1.iterrows():
            if continue_outer_loop:  # Reset the flag's value before the inner loop
                continue_outer_loop = False
            if row1_single['matched']:  # Skip rows that have already been matched
                continue
            chi1_1, chi2_1 = row1_single['chi1'], row1_single['chi2']
        
            for index2, row2_single in df2.iterrows():
                if row2_single['matched']:  # Skip rows that have already been matched
                    continue
                chi1_2, chi2_2 = row2_single['chi1'], row2_single['chi2']

                # Check if the angles of row1 and row2 match
                
                if row1_single['rotamer'][:2] == row2_single['rotamer'][:2]: #angle_difference_greater_than_15(chi1_1, chi1_2):
                        match_count = 1
                        df1.loc[index1, 'matched'] = True
                        df2.loc[index2, 'matched'] = True
                        continue_outer_loop = True
                        break
                        # If there's a match and more than one row1_single, move on to the next row1_single
                else:
                        no_match = 1
        
        # Save the evaluation results
            if no_match == 0:
                TP_count = 1
            elif match_count == 0:
                FP_count = 1
            elif no_match > 0 and match_count > 0:
                if len(df1) > len(df2):
                    loss_count = 1
                elif len(df2) > len(df1):
                    gain_count = 1
                else:
                    share = 1


            # Return the evaluation results
            return {'TP': TP_count, 'loss': loss_count, 'gain': gain_count, 'FP': FP_count, 'share_diff': share}


def calc_sum(results_df):
    # Calculate the sum of TP, FP, Gain, and Loss columns
    total_count = results_df["TP"].sum() + results_df["FP"].sum() + results_df["gain"].sum() + results_df["loss"].sum() + results_df["share_diff"].sum()
    # Calculate the percentages of TP, FP, Gain, and Loss
    results_df["TP_pct"] = results_df["TP"] / total_count * 100
    results_df["FP_pct"] = results_df["FP"] / total_count * 100
    results_df["gain_pct"] = results_df["gain"] / total_count * 100
    results_df["loss_pct"] = results_df["loss"] / total_count * 100
    results_df["share_diff"] = results_df["share_diff"] / total_count * 100

    print("Total count:", total_count)
    print("TP count:", results_df["TP"].sum())
    print("FP count:", results_df["FP"].sum())
    print("Gain count:", results_df["gain"].sum())
    print("Loss count:", results_df["loss"].sum())
    print("Share Some count:", results_df["share_diff"].sum())
    

    # Print the percentages
    print("TP %:", results_df["TP_pct"].sum())
    print("FP %:", results_df["FP_pct"].sum())
    print("Gain %:", results_df["gain_pct"].sum())
    print("Loss %:", results_df["loss_pct"].sum())
    print("Share Some %:", results_df["share_diff"].sum())

    all_array = [results_df["TP_pct"].sum(), results_df["FP_pct"].sum(), results_df["gain_pct"].sum(),results_df["loss_pct"].sum(), results_df["share_diff"].sum()]
    return all_array


def plot_results(all_, alt_):
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
        # Set y-ticks and x-labels font size to 20
    ax1.tick_params(axis='y', labelsize=14)
    ax2.tick_params(axis='x', labelsize=14)
    ax1.tick_params(axis='x', labelsize=14)
    # Set x-labels font size to 20
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[0], color='#191970', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[1], bottom=alt_[0], color='#097969', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[2], bottom=alt_[0] + alt_[1], color='#FF5733', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[3], bottom=alt_[0] + alt_[1] + alt_[2], color='#C70039', width=bar_width)
    alt_bar = ax2.bar('Deposited Alternative Conformer', alt_[4], bottom=alt_[0] + alt_[1] + alt_[2] + alt_[3], color='#FFC300', width=bar_width)

    # Adjust the space between the two bar plots
    plt.subplots_adjust(wspace=space_between)

    # Add a common y-axis label
    fig.text(0.02, 0.5, 'Percentage', va='center', rotation='vertical', fontsize=16)

    # Save the bar plot to a file
    plt.savefig(f'RotamerStatus_stacked_bars.png')
    
def save_legend():
    print('save legend')
    # Define the colors and labels for the legend
    colors = ['#191970', '#097969', '#FF5733', '#C70039', '#FFC300']
    labels = ['Consistent Rotamer(s)', 'Different Rotamer(s)', 'Additional Rotamer(s) in qFit', 'Additional Rotamer(s) in Deposited', 'Consistent & Different Rotamers']

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Create the legend handles with colored boxes
    handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors]

    # Create the legend with the handles and labels
    legend = ax.legend(handles, labels, loc='upper right')

    # Set the title for the legend
    legend.set_title('Legend')
    fig.savefig('/Users/stephaniewanko/Downloads/temp/qfit_test_set/RotamerStatus_legend.png', bbox_inches='tight', dpi=600)


def main():
    qfit = read_files('/Users/stephaniewanko/Downloads/temp/qfit_test_set', 'final2', 58,62)
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
        #unique_residues = [5,27]
        for resi in unique_residues:
            # Filter the dataframes to only include rows with the current residue
            df1_residue = df1_filtered[df1_filtered['resi'] == resi]
            df2_residue = df2_filtered[df2_filtered['resi'] == resi]

            # Subset df1_residue and df2_residue to chain, resi, PDB, rotamer, chi1, chi2
            df1_residue = df1_residue[['chain', 'resi', 'PDB', 'rotamer', 'chi1', 'chi2']]
            df2_residue = df2_residue[['chain', 'resi', 'PDB', 'rotamer', 'chi1', 'chi2']]

                        # Remove duplicate rows from df1_residue and df2_residue if rotamers are the same
            df1_residue = df1_residue.drop_duplicates(subset=['rotamer'])
            df2_residue = df2_residue.drop_duplicates(subset=['rotamer'])


            # Compare the rotamers of the current residue in both dataframes
            comparison_result = compare_rotamers(df1_residue, df2_residue)
            comparison_result['residue'] = resi
            comparison_result['PDB'] = pdb

            # Check if df1_residue has more than 1 line
            if len(df1_residue) > 1:
                results.append(comparison_result)
                results_alt.append(comparison_result)
            else:
                results.append(comparison_result)

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
    print('All:')
    all_ = calc_sum(results_df)
    print('Alt:')
    alt_ = calc_sum(results_df_alt)  

    save_legend()

    # Save the DataFrame to a CSV file
    results_df.to_csv("comparison_results.csv", index=False)
    results_df_alt.to_csv("comparison_results_alt.csv", index=False)


    # Call the plot_results function
    plot_results(all_, alt_)


        

if __name__ == "__main__":
    main()
