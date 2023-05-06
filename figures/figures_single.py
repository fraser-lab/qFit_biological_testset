import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import glob
import os
import matplotlib.patches as mpatches
import argparse



pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)

def parse_args():
    parser = argparse.ArgumentParser(description='Compare the results of different water picking methods')
    parser.add_argument('--base_dir', type=str, help='The base directory of the water picking results')
    parser.add_argument('--comparison', type=str, help='The name of the comparison method')
    return parser.parse_args()

def rotamer_subset_compare(subset_rotamer, comparison):
    rotamer_ = []
    comparison_rotamer = comparison + '_rotamer'
    for i in subset_rotamer['original_PDB'].unique():
        tmp = subset_rotamer[subset_rotamer['original_PDB'] == i]
        for r in tmp['chain_resi'].unique():
            num_original = len(set(tmp[tmp['chain_resi'] == r][f'original_rotamer']))
            num_compairson = len(set(tmp[tmp['chain_resi'] == r][comparison_rotamer]))
            if set(tmp[tmp['chain_resi']==r]['original_rotamer']) - set(tmp[tmp['chain_resi']==r][comparison_rotamer]):
                rotamer = 'Different'
            #elif bool(set(tmp[tmp['chain_resi']==r]['original_rotamer']) & set(tmp[tmp['chain_resi']==r][comparison_rotamer])) == True:
            if len(set(tmp[tmp['chain_resi']==r]['original_rotamer'])) > len(set(tmp[tmp['chain_resi']==r][comparison_rotamer])):
                    rotamer = 'Gain in Original'
            elif len(set(tmp[tmp['chain_resi']==r]['original_rotamer'])) < len(set(tmp[tmp['chain_resi']==r][comparison_rotamer])):
                    rotamer = 'Gain in Comparison'
            else:
                rotamer = 'Same'
            rotamer_.append(tuple((i, r, rotamer, num_original, num_compairson)))
    rotamer_comp = pd.DataFrame(rotamer_, columns =['PDB', 'chain_resi', 'Rotamer', 'Num_Original', 'Num_Comp'])
    return rotamer_comp


def read_files(base_dir, comparison, file_start, file_end):
    #bfactors
    all_files = glob.glob(base_dir + '/' + comparison + "/*B_factors.csv")
    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, sep=',', header=0)
        df['PDB'] = filename[file_start:file_end] #file_start:file_end
        li.append(df)
    b_factor = pd.concat(li, axis=0, ignore_index=True)
    b_factor['category'] = comparison

    #rotamers
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

    #rvalues
    all_files = glob.glob(base_dir + '/' + comparison +  "/*_rvalues.csv")
    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col = None, sep=',', header=0)
        df['PDB'] = filename[file_start:file_end]
        li.append(df)
    r_values = pd.concat(li, axis=0, ignore_index=True)
    return b_factor, r_values, rotamer


def compare_rvalues(comparison_rvalues, original_rvalues, comparison, base_dir):
    original_rvalues.columns = 'original_' + original_rvalues.columns.values
    comparison_rvalues.columns = 'b_factor' + '_' + comparison_rvalues.columns.values

    comparison_Rfree = 'b_factor' + '_Rfree_qFit'
    r_values = comparison_rvalues #original_rvalues.merge(comparison_rvalues, left_on=['original_PDB'], right_on=[comparison_PDB])
    #calculating the difference in rvalues between each step
    r_values['Rfree_Differences'] = r_values['b_factor_Rfree_pre'] - r_values[comparison_Rfree]
    r_values['Rgap_deposited'] = r_values['b_factor_Rwork_pre'] - r_values['b_factor_Rfree_pre']
    r_values['Rgap_qFit'] = r_values['b_factor_Rwork_qFit'] - r_values['b_factor_Rfree_qFit']
    print('rvalues outliers:')
    print(r_values[r_values['b_factor_Rfree_qFit']> 0.25])
    r_values['Rfree_improve'] = np.where(((r_values['b_factor_Rfree_pre'] - r_values[comparison_Rfree]) > 0), 1, 0)
    print((len(r_values[r_values['Rfree_improve']==1].index))/len(r_values.index))
    print(r_values['Rfree_Differences'].mean())


    #labeling structures that will be removed due to high rfree values
    r_values['Rfree_YN'] = np.where(((r_values['b_factor_Rfree_pre'] - r_values[comparison_Rfree]) < -0.025), 1, 0)
    #plotting the difference in rvalues
        # Calculate the range of Rfree values for b_factor_Rfree_pre and b_factor_Rfree_post
    rfree_pre_min = r_values['b_factor_Rfree_pre'].min()
    rfree_pre_max = r_values['b_factor_Rfree_pre'].max()
    rfree_post_min = r_values[comparison_Rfree].min()
    rfree_post_max = r_values[comparison_Rfree].max()

    print(f"Min Rfree value for b_factor_Rfree_pre: {rfree_pre_min}")
    print(f"Max Rfree value for b_factor_Rfree_pre: {rfree_pre_max}")
    print(f"Min Rfree value for b_factor_Rfree_post: {rfree_post_min}")
    print(f"Max Rfree value for b_factor_Rfree_post: {rfree_post_max}")
    
    fig=plt.figure()
    sns.lmplot('b_factor_Rfree_pre', comparison_Rfree, data=r_values, fit_reg=False)
    slope = 1
    y_intercept = 0
    x_fit = np.linspace(0.12, 0.27, 100)
    y_fit = slope * x_fit + y_intercept
    plt.plot(x_fit, y_fit, color='red')
    plt.xlabel('Re-Refined Deposited Rfree', fontsize=12)
    plt.text(0.12, 0.23, 'Better Deposited Rfree', fontsize=12) 
    plt.text(0.19, 0.15, f'Better qFit Rfree', fontsize=12)
    plt.ylabel('qFit Rfree', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0.12, 0.27) # Set the same range for x-axis
    plt.ylim(0.12, 0.27) # Set the same range for y-axis
    plt.savefig(base_dir + '/Refinement_scatterplot_single.png')


    fig = plt.figure()
    sns.kdeplot(r_values['Rfree_Differences'], shade=True)
    plt.axvline(x = 0, color = 'red')
    plt.xlabel('Rfree Difference', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(base_dir + f'/Rfree_Differences_distribution_single.png')

    fig = plt.figure()
    sns.kdeplot(r_values['Rgap_deposited'], shade=True, color='#FF5733', label='Deposited')
    sns.kdeplot(r_values['Rgap_qFit'], shade=True, color='#097969', label='qFit')
    plt.xlabel('R-Gap (Rwork-Rfree)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.savefig(base_dir + f'/Rgap_distribution_single.png')

    print(f'Deposited Median:')
    print(r_values['Rgap_deposited'].median())
    print(f'qFit Median:')
    print(r_values['Rgap_qFit'].median())

    print(r_values['Rfree_Differences'].median())
    print(len(r_values[r_values['Rfree_Differences']>=0].index))


def compare_rotamer(original_rotamer, comparison_rotamer, comparison, base_dir):    
    original_rotamer.columns = 'original_' + original_rotamer.columns.values
    comparison_rotamer.columns = comparison + '_' + comparison_rotamer.columns.values

    comparison_PDB = comparison + '_PDB'
    comparison_chain = comparison + '_chain'
    comparison_resi = comparison + '_resi'
    rotamer_all = original_rotamer.merge(comparison_rotamer, left_on=['original_PDB', 'original_chain', 'original_resi'], right_on=[comparison_PDB, comparison_chain, comparison_resi])
    rotamer_all['chain_resi'] = rotamer_all['original_chain'] + '_' + rotamer_all['original_resi'].astype(str)
    comp_rotamer= rotamer_subset_compare(rotamer_all, comparison)

    comp_rotamer.to_csv(f'{base_dir}/comp_rotamer.csv')

    comp_eval = comparison + '_evaluation'
    sns.countplot(x=comparison_PDB, data=comparison_rotamer[comparison_rotamer[comp_eval] == "OUTLIER"])
    plt.ylabel('PDB Id')
    plt.xlabel('Number of OUTLIER Rotamers')
    plt.savefig(base_dir + f'/Outlier_Rotamers_{comparison}.png')
    comp_rotamer = comp_rotamer[comp_rotamer['Num_Comp']>1]

    different = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Different'].index)
    same = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Same'].index)
    gain_comp = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Gain in Comparison'].index)
    gain_original = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Gain in Original'].index)
    print(f'% Different:{different/len(comp_rotamer.index)}')
    print(f'% Same:{same/len(comp_rotamer.index)}')
    print(f'% Gain {comparison}:{gain_comp/len(comp_rotamer.index)}')
    print(f'% Gain Original: {gain_original/len(comp_rotamer.index)}')
    print(f'Gain Original: {gain_original}')
    print(f'All: {len(comp_rotamer.index)}')
    all_ = [same, gain_comp, gain_original]
    labels_ = ['Same', f'Additional rotamer unique to qFit model', f'Additional rotamer unique to Deposited model']

    fig = plt.figure()
    ax = plt.bar([1,2,3], all_, color = ['#191970', '#097969', '#FF5733'])
    plt.tight_layout()
    plt.xticks([1,2,3], labels_, rotation = 10) #, color=['black', 'red', 'blue', 'cyan']
    plt.ylabel('Number of Residues', fontsize=12)
    plt.savefig(f'{base_dir}/RotamerStatus_{comparison}.png', bbox_inches="tight")

    #compare rotamer of only those with altlocs in deposited PDB
    comp_rotamer_alt = comp_rotamer[comp_rotamer['Num_Original']>1]
    comp_rotamer_alt.to_csv(f'{base_dir}/{comparison}_comp_rotamer_alt.csv')


    print('num residues with altloc in qfit')
    print(len(comp_rotamer[comp_rotamer['Num_Comp']>1].index))
    print(f'num residue with alt loc {len(comp_rotamer_alt.index)}')

    different = len(comp_rotamer_alt[comp_rotamer_alt['Rotamer'] == 'Different'].index)
    same = len(comp_rotamer_alt[comp_rotamer_alt['Rotamer'] == 'Same'].index)
    gain_comp = len(comp_rotamer_alt[comp_rotamer_alt['Rotamer'] == 'Gain in Comparison'].index)
    gain_original = len(comp_rotamer_alt[comp_rotamer_alt['Rotamer'] == 'Gain in Original'].index)
    print(f'% Different:{different/len(comp_rotamer_alt.index)}')
    print(f'% Same:{same/len(comp_rotamer_alt.index)}')
    print(f'% Gain {comparison}:{gain_comp/len(comp_rotamer_alt.index)}')
    print(f'% Gain Original: {gain_original/len(comp_rotamer_alt.index)}')
    print(f'Gain Original: {gain_original}')
    all_alt = [same, gain_comp, gain_original]
    array1 = np.array(all_)
    array2 = np.array(all_alt)

    # Calculate the percentage of total value for each array
    array1_percentage = array1 / np.sum(array1) * 100
    array2_percentage = array2 / np.sum(array2) * 100

    fig = plt.figure()
    # Create the stacked bar plot of array1_percentage with each value being the bottom of the next

    # Create a numpy array for each set of values, and transpose them
    all_array = np.array([all_])
    altlocs_array = np.array([all_alt])
    array1_percentage = all_array / np.sum(array1) * 100
    array2_percentage = altlocs_array / np.sum(array2) * 100
    all_array = array1_percentage.T
    altlocs_array = array2_percentage.T


    # Create    # Create a figure with two subplots, side by side
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)

    # Create a stacked bar plot for the first set of values
    ax1.bar('All', all_array[0], color='#191970', width=0.5)
    ax1.bar('All', all_array[1], bottom=all_array[0], color='#097969', width=0.5)
    ax1.bar('All', all_array[2], bottom=all_array[0] + all_array[1], color='#FF5733', width=0.5)

    # Create a stacked bar plot for the second set of values
    ax2.bar('Deposited AltLocs', altlocs_array[0], color='#191970', width=0.5)
    ax2.bar('Deposited AltLocs', altlocs_array[1], bottom=altlocs_array[0], color='#097969', width=0.5)
    ax2.bar('Deposited AltLocs', altlocs_array[2], bottom=altlocs_array[0] + altlocs_array[1], color='#FF5733', width=0.5)

    # Add a common y-axis label
    fig.text(0.04, 0.5, 'Percentage', va='center', rotation='vertical', fontsize=12)

    plt.savefig(f'{base_dir}/RotamerStatus_stacked_half_width.png')

    fig = plt.figure()
    colors = ['#191970', '#097969', '#FF5733']
    f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
    handles = [f("s", colors[i]) for i in range(3)]
    labels = ['Same', 'Additional rotamer unique to qFit model', 'Additional rotamer unique to Deposited model']
    legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False)

    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(f'{base_dir}/RotamerStatus_legend.png', dpi="figure", bbox_inches=bbox)


def compare_bfactor(original_bfactor, comparison_bfactor, comparison, base_dir): 
    original_bfactor.columns = 'original_' + original_bfactor.columns.values
    comparison_bfactor.columns = 'b_factor' + '_' + comparison_bfactor.columns.values
    comparison = 'b_factor'
    comparison_PDB = comparison + '_PDB'
    comparison_chain = comparison + '_chain'
    comparison_resi = comparison + '_resi'
    comparison_b = 'b_factor' + '_b_factor'
    comparison_altloc = 'b_factor' + '_num_altlocs'

    print('number of residues:')
    print(len(original_bfactor.index))
    print('number alt loc in original:')
    print(len(original_bfactor[original_bfactor['original_num_altoc']>1].index))
    print('number alt loc in qFit:')
    print(comparison_bfactor.head())
    print(len(comparison_bfactor[comparison_bfactor[comparison_altloc]>1].index))
    b_factor_all = original_bfactor.merge(comparison_bfactor, left_on=['original_PDB', 'original_chain', 'original_resi'], right_on=[comparison_PDB, comparison_chain, comparison_resi])
    b_factor_all['bfactor_diff'] =  b_factor_all[comparison_b] - b_factor_all['original_b_factor']
    b_factor_all['altloc_diff'] = b_factor_all[comparison_altloc]- b_factor_all['original_num_altoc']
    print('b-factors')
    print(len(np.unique(b_factor_all['original_PDB'])))
    b_factor_all =b_factor_all[b_factor_all['altloc_diff']<5]
    #compare altlocs
    fig = plt.figure()
    ax = sns.countplot(x='altloc_diff', data=b_factor_all, color='#24248f')
    #ax = sns.histplot(b_factor_all['altloc_diff'], bins=[-4, -3, -2, -1, 0, 1, 2, 3, 4], shrink=.8, color='#24248f')
    for p in ax.patches:
        height = p.get_height() # get the height of each bar
        # adding text to each bar
        ax.text(x = p.get_x()+(p.get_width()/2), y = height+0.3, s = '{:.0f}'.format(height), ha = 'center') 
    #plt.xticks(arr_div, arr_div_r, fontsize=12)
    plt.xlabel('Difference in Number of Alternative Conformers', fontsize=12)
    plt.ylabel('Number of Residues', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    #plt.text(-4.8, 6000, f'More Altlocs in qFit Model') 
    #plt.text(2, 6000, 'More Altlocs in Deposited Model') 
    plt.savefig(f'{base_dir}/altloc_diff.png', bbox_inches="tight")

    fig = plt.figure()
    sns.boxplot(b_factor_all['altloc_diff'])
    #plt.text(-4.8, 6000, f'More Alt locs with {comparison}') 
    #plt.text(2, 6000, 'More Altlocs with Original') 
    plt.savefig(f'{base_dir}/altloc_diff_box.png')

    #compare bfactor
    fig = plt.figure()
    #remove zeros 
    #b_factor_zeros = b_factor_all[b_factor_all['bfactor_diff']!=0]
    #print(b_factor_zeros['bfactor_diff'].median())
    #b_factor_zeros.to_csv('b_factor_zeros.csv')
    sns.histplot(b_factor_all['bfactor_diff'], bins=30)
    plt.axvline(x = 0, color = 'r', label = '')
    #plt.text(-250, 0.01, f'Higher B-factor in qFit Model', fontsize=12) 
    #plt.text(0, 0.01, f'Higher B-factors in Deposited Model', fontsize=12)
    plt.xlabel('Difference in Residue B-factors', fontsize=12)
    plt.ylabel('Number of Residues', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(f'{base_dir}/bfactor_diff.png')

    print(f'Average B-factor difference: {b_factor_all["bfactor_diff"].mean()}')
    print(f'Median Altloc difference: {b_factor_all["altloc_diff"].median()}')
    print(f'Mean Altloc difference: {b_factor_all["altloc_diff"].mean()}')

    # print('The following residues gain more than two altloc, you should look at them:')
    # print(b_factor_all[b_factor_all['altloc_diff'] > 2][['original_PDB', 'original_chain', 'original_resi', 'original_num_altoc', comparison_altloc]])

    # print('The following residues lose more than two altloc, you should look at them:')
    # print(b_factor_all[b_factor_all['altloc_diff'] < -2][['original_PDB', 'original_chain', 'original_resi', 'original_num_altoc', comparison_altloc]])


    # print('The following residues have a bfactor difference greater than 100, you should look at them:')
    # print(b_factor_all[b_factor_all['bfactor_diff'] > 100][['original_PDB', 'original_chain', 'original_resi', 'original_b_factor', comparison_b]])
    # print(b_factor_all[b_factor_all['bfactor_diff'] < -100][['original_PDB', 'original_chain', 'original_resi', 'original_b_factor', comparison_b]])


def main():
    args = parse_args()
    #read in original
    original_bfactor, original_rvalues, original_rotamer = read_files(args.base_dir, 'single_conf', 63, 67)  #51, 55
    #read in comparison
    comparison_bfactor, comparison_rvalues, comparison_rotamer = read_files(args.base_dir, args.comparison, 57,61) #read in cctbx: 65-69 b-factor:60, 64 BIC: 55,59 BIC b: 57,61
    print(comparison_rotamer.head())
    compare_rvalues(comparison_rvalues, original_rvalues, args.comparison, args.base_dir)
    compare_rotamer(original_rotamer, comparison_rotamer, args.comparison, args.base_dir)
    compare_bfactor(original_bfactor, comparison_bfactor, args.comparison, args.base_dir)
if __name__ == '__main__':
    main()
