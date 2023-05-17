import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns


def parse_molprobity_scores(file_path):
    with open(file_path, 'r') as molfile:
        towrite = ""
        for line in molfile:
            if line.startswith("  Ramachandran outliers ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                towrite+=lineparts[1][0:-1].strip()
            elif line.startswith("                favored ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                towrite+=lineparts[1][0:-1].strip()
            elif line.startswith("  Rotamer outliers      ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                towrite+=lineparts[1][0:-1].strip()
            elif line.startswith("  C-beta deviations     ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                towrite+=lineparts[1].strip()
            elif line.startswith("  Clashscore            ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                if "percentile" in line:
                    score=lineparts[1].split("(")[0].strip()
                    towrite+=score
                else:
                    towrite+=lineparts[1].strip()
            elif line.startswith("  RMS(bonds)            ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                towrite+=lineparts[1].strip()
            elif line.startswith("  RMS(angles)           ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                towrite+=lineparts[1].strip()
            elif line.startswith("  MolProbity score      ="):
                lineparts=line[0:-1].split("=")
                towrite+=","
                if "percentile" in line:
                    score=lineparts[1].split("(")[0].strip()
                    towrite+=score
                else:
                    towrite+=lineparts[1].strip()
    return towrite




def read_files(directory_path):
    files = [f for f in os.listdir(directory_path) if f.endswith("validation.txt")]
    data = []
    for file in files:
        file_path = os.path.join(directory_path, file)
        molprobity_scores = parse_molprobity_scores(file_path)
        PDB = file[:4]
        molprobity_scores = PDB + "," + molprobity_scores
        data.append(molprobity_scores.split(","))

    df = pd.DataFrame(data, columns=["PDB", "", "Ramachandran outliers", "Ramachandran favored", "Rotamer outliers", "C-beta deviations", "Clashscore", "RMS(bonds)", "RMS(angles)", "MolProbity score"])

    # loop through each column in df[2:] and make the column type numeric
    df.iloc[:, 2:] = df.iloc[:, 2:].apply(pd.to_numeric, errors='coerce')
    return df


directory_path_multi = "/Users/stephaniewanko/Downloads/temp/qfit_test_set/final2/"
directory_path_single = "/Users/stephaniewanko/Downloads/temp/qfit_test_set/single_conf/"
df = read_files(directory_path_multi)
single_molprobity = read_files(directory_path_single)




# Use sns.histplot instead of sns.kdeplot for histogram
for col in df.columns[2:]:
    df[col] = pd.to_numeric(df[col], errors='coerce')
    single_molprobity[col] = pd.to_numeric(single_molprobity[col], errors='coerce')
    # create a histplot of the column
    # set the title of the histogram to the name of the column
        # Determine the range of values for both datasets
    min_value = min(df[col].min(), single_molprobity[col].min())
    max_value = max(df[col].max(), single_molprobity[col].max())

    # Create a common bin size for both histograms
    bins = np.linspace(min_value, max_value, 13)


    hist_single = sns.histplot(x=col, data=single_molprobity, color='darkgreen', bins=bins, kde=False, alpha=1.0, label='qFit Model')
    hist_df = sns.histplot(x=col, data=df, color='darkmagenta', kde=False, bins=bins, alpha=0.8, label='Deposited Model')
    plt.legend()
    # set x-axis label
    plt.xlabel(col)
    # set y-axis label
    plt.ylabel("Count")
    # Save the stacked histplot as a figure with the name of the column + stacked_histplot
    plt.savefig(col + "_stacked_histplot.png")
    # Clear the current figure
    plt.clf()
    
# Calculate the median, Q1, and Q3 for each column in the qFit and single df
qfit_summary = df.iloc[:, 2:].describe().loc[['25%', '50%', '75%']]
single_summary = single_molprobity.iloc[:, 2:].describe().loc[['25%', '50%', '75%']]

# Rename the index for better readability
qfit_summary.index = ['Q1 (qFit)', 'Median (qFit)', 'Q3 (qFit)']
single_summary.index = ['Q1 (Single)', 'Median (Single)', 'Q3 (Single)']

# Combine the two summary dataframes
combined_summary = pd.concat([qfit_summary, single_summary])

# Print the combined summary dataframe
print(combined_summary)
