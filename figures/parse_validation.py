import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind


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



directory_path_new = "/final2/"
directory_path_base = "/single_conf/"
new = read_files(directory_path_new)
base = read_files(directory_path_base)
base.to_csv('base_qFit_validation.csv')
new.to_csv('new_qFit_validation.csv')


plt.rcParams.update({'font.size': 16})
sns.set_style("white")
angstrom_symbol = "\u00C5"
sigma_symbol = "\u03A3"
beta_symbol = "\u03B2"
degree_symbol = '\u00b0'
plt.rcParams["figure.dpi"] = 500

for col in df.columns[2:]:
    new[col] = pd.to_numeric(new[col], errors='coerce')
    base[col] = pd.to_numeric(base[col], errors='coerce')
    
    # Determine the range of values for both datasets
    min_value = min(new[col].min(), base[col].min())
    max_value = max(new[col].max(), base[col].max())
    

    # Run a ttest between the values in the df[col] and the single_molprobity[col]
    t_stat, p_val = ttest_ind(new[col].dropna(), base[col].dropna())
    print(f'T-test results for {col}: T-Stat = {t_stat}, P-value = {p_val}')
    

    # Create a common bin size for both histograms
    p1 = plt.figure()
    bins = np.linspace(min_value, max_value, 13)
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharey=True)  # sharex=True, Add sharey=True to make the y-axis the same height


    # Flip the order of the histograms
    ax1, ax2 = ax2, ax1  
    # Remove the border between the graphs
    fig.subplots_adjust(hspace=0)
    ax1.set_xlim(min_value, max_value)
    ax2.set_xlim(min_value, max_value)
    ax2.set_xticks([])  # Remove x-axis ticks for the top subplot
    
    sns.histplot(x=col, data=df, color='darkmagenta', kde=False, bins=bins, alpha=1.0, ax=ax1)
    ax1.set_ylabel("# of Structures")
    xlab = f'{col}'
    ax1.legend(["New qFit Model"])
    ax1.set_xlabel(xlab)
    #ax1.set_xticklabels([])  # Remove x-axis labels for the top subplot

    # Create a histplot of the column for Deposited Model
    sns.histplot(x=col, data=single_molprobity, color='darkgreen', kde=False, bins=bins, alpha=1.0, ax=ax2)
    ax2.set_ylabel("# of Structures")

    ax2.legend(["Baseline qFit Model"])

    # Save the stacked histplot as a figure with the name of the column + stacked_histplot
    # Increase the resolution of the figures
    # Ensure the x-axis is not cut off
    plt.tight_layout()
    
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
combined_summary.to_csv('combined_summary.csv')
