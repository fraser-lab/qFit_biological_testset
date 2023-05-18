#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:47:00 2023

@author: ashraya
"""

import numpy as np
from qfit import Structure
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


base_dir="/Users/ashraya/qfit_summit/summit_2/macdomain_maps/"

def plot_qscore():
    qscore_df=pd.read_csv(base_dir+"correct_qscore_fraction.csv")
    ax1=sns.lineplot(data=qscore_df, x="resolution", y="correct qscore",ci='sd')
    plt.xlim(0.8,3.0) 
    ax1.set_xticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0])           
    plt.xlabel("Resolution ($\AA$)",fontsize=18)
    plt.ylabel("Fraction of residues matching\nin Q-score with ground truth",fontsize=18)       
    plt.savefig(base_dir+"qscore_comparison_plot.png",dpi=300,bbox_inches="tight")   
        

number_of_runs=10
r10s=range(8,31)
opfile=open(base_dir+"correct_qscore_fraction.csv","w")

opfile.write("resolution,run,correct qscore\n")
for r10 in r10s:
    resol = r10 * 0.1
    resol = str(round(resol,1))
    print(resol)
    for run in range(1,number_of_runs+1):
        correct_qscore=0
        max_qscore_dict=dict()
        ref_struct=Structure.fromfile(base_dir+"xray_noise_modeling/p1_box_pdbs/7kr0_shift_bfac_"+resol+"_p1.pdb__Q__7kr0_resol_"+resol+"_shake_noisy_run_"+str(run)+"_2mFo-DFc.ccp4.pdb").reorder()
        for ref_residue in (
                ref_struct.extract("record", "ATOM")  # Don't analyse metals/ligands
                .extract("resn", "HOH", "!=")  # Don't analyse waters
                .extract(
                    "name", "H", "!="
                )  # Sometimes backbone N-H atoms are present in some altlocs, not all. Avoid analysing them.
                .extract(
                    "e", "H", "!="
                )  # Sometimes His protonation states differ between altlocs. Avoid analysing all H.
        ).residue_groups:
            altlocs = sorted(list(set(ref_residue.altloc)))
            max_qscore=-999.0
            if '' in altlocs and len(altlocs) > 1:
                read_sidechain=False
                main_chain_qs=[]
                for altloc in altlocs:
                    alt_qscores=[]
                    ref_conf=ref_residue.extract("altloc",altloc)
                    if altloc=='':
                        for b in ref_conf.b:
                            main_chain_qs.append(b)      
                    else:
                        for b in ref_conf.b:
                            alt_qscores.append(b)
                        read_sidechain=True
                    if read_sidechain:
                        alt_avg_q=(sum(main_chain_qs)+sum(alt_qscores))/float(len(main_chain_qs)+len(alt_qscores))
                        if alt_avg_q > max_qscore:
                            max_qscore=alt_avg_q
            else:
                 for altloc in altlocs:
                     alt_qscores=[]
                     ref_conf=ref_residue.extract("altloc",altloc)
                     for b in ref_conf.b: 
                         alt_qscores.append(b)
                     alt_avg_q=np.mean(alt_qscores)
                     if alt_avg_q > max_qscore:
                            max_qscore=alt_avg_q
            max_qscore_dict[ref_residue.resi[0]]=max_qscore
        qfit_struct=Structure.fromfile(base_dir+"resol_"+resol+"/run_"+str(run)+"_final/qfit_bic/7kr0_resol_"+resol+"_shake_noisy_run_"+str(run)+"_qFit.pdb__Q__7kr0_resol_"+resol+"_shake_noisy_run_"+str(run)+"_2mFo-DFc.ccp4.pdb").reorder()
        for ref_residue in (
                qfit_struct.extract("record", "ATOM")  # Don't analyse metals/ligands
                .extract("resn", "HOH", "!=")  # Don't analyse waters
                .extract(
                    "name", "H", "!="
                )  # Sometimes backbone N-H atoms are present in some altlocs, not all. Avoid analysing them.
                .extract(
                    "e", "H", "!="
                )  # Sometimes His protonation states differ between altlocs. Avoid analysing all H.
        ).residue_groups:
            resi=ref_residue.resi[0]

            altlocs = sorted(list(set(ref_residue.altloc)))
            max_qscore=-999.0
            if '' in altlocs and len(altlocs) > 1:
                read_sidechain=False
                main_chain_qs=[]
                for altloc in altlocs:
                    alt_qscores=[]
                    ref_conf=ref_residue.extract("altloc",altloc)
                    if altloc=='':
                        for b in ref_conf.b:
                            main_chain_qs.append(b)      
                    else:
                        for b in ref_conf.b:
                            alt_qscores.append(b)
                        read_sidechain=True
                    if read_sidechain:
                        alt_avg_q=(sum(main_chain_qs)+sum(alt_qscores))/float(len(main_chain_qs)+len(alt_qscores))
                        if alt_avg_q > max_qscore:
                            max_qscore=alt_avg_q
            else:
                 for altloc in altlocs:
                     alt_qscores=[]
                     ref_conf=ref_residue.extract("altloc",altloc)
                     for b in ref_conf.b: 
                         alt_qscores.append(b)
                     alt_avg_q=np.mean(alt_qscores)
                     if alt_avg_q > max_qscore:
                            max_qscore=alt_avg_q

            
            qscore_diff=max_qscore-max_qscore_dict[resi]

            if qscore_diff >= -0.01:
                correct_qscore+=1
        total_correct_frac=correct_qscore/169.0

        opfile.write(resol+","+str(run)+","+str(round(total_correct_frac,2))+"\n")
opfile.close()
plot_qscore()
