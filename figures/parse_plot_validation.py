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


directory_path_multi = "/Users/stephaniewanko/Downloads/temp/qfit_test_set/final/"
directory_path_single = "/Users/stephaniewanko/Downloads/temp/qfit_test_set/single_conf/"
df = read_files(directory_path_multi)
single_molprobity = read_files(directory_path_single)




# Use sns.histplot instead of sns.kdeplot for histogram
for col in df.columns[2:]:
    df[col] = pd.to_numeric(df[col], errors='coerce')
    single_molprobity[col] = pd.to_numeric(single_molprobilty[col], errors='coerce')
    # create a histplot of the column
    hist_df = sns.histplot(x=col, data=df, color='blue', kde=False, bins=12, alpha=0.5)
    hist_single = sns.histplot(x=col, data=single_molprobity, color='green', kde=False, bins=12, alpha=0.5)

    # set x-axis label
    plt.xlabel(col)
    # set y-axis label
    plt.ylabel("Count")

    # Add legend
    plt.legend(["qFit Model", "Deposited Model"])

    #  save the histplot as a figure with the name of the column + histplot
    plt.savefig(col + "_histplot.png")
    # clear the current figure
    plt.clf()
    
