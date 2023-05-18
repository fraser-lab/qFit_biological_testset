#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 14:19:16 2023

@author: ashraya
"""
from qfit import Structure
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

base_dir="/Users/ashraya/qfit_summit/summit_2/macdomain_maps/"

def my_rmsd(res1, res2):
        coor1 = res1.coor
        coor2 = res2.coor

        if coor1.shape != coor2.shape:
            raise ValueError("Coordinate shapes are not equivalent")
        if "TYR" in res1.resn:
            idx_cd1 = res2.name.tolist().index("CD1")
            idx_cd2 = res2.name.tolist().index("CD2")
            idx_ce1 = res2.name.tolist().index("CE1")
            idx_ce2 = res2.name.tolist().index("CE2")
            coor3 = np.copy(coor2)
            coor3[idx_cd1], coor3[idx_cd2] = coor2[idx_cd2], coor2[idx_cd1]
            coor3[idx_ce1], coor3[idx_ce2] = coor2[idx_ce2], coor2[idx_ce1]
            diff = (coor1 - coor2).ravel()
            diff2 = (coor1 - coor3).ravel()
            return min(
                np.sqrt(3 * np.inner(diff, diff) / diff.size),
                np.sqrt(3 * np.inner(diff2, diff2) / diff2.size),
            )
        if "PHE" in res1.resn:
            idx_cd1 = res2.name.tolist().index("CD1")
            idx_cd2 = res2.name.tolist().index("CD2")
            idx_ce1 = res2.name.tolist().index("CE1")
            idx_ce2 = res2.name.tolist().index("CE2")
            coor3 = np.copy(coor2)
            coor3[idx_cd1], coor3[idx_cd2] = coor2[idx_cd2], coor2[idx_cd1]
            coor3[idx_ce1], coor3[idx_ce2] = coor2[idx_ce2], coor2[idx_ce1]
            diff = (coor1 - coor2).ravel()
            diff2 = (coor1 - coor3).ravel()
            return min(
                np.sqrt(3 * np.inner(diff, diff) / diff.size),
                np.sqrt(3 * np.inner(diff2, diff2) / diff2.size),
            )
        if "ARG" in res1.resn:
            idx_nh1 = res2.name.tolist().index("NH1")
            idx_nh2 = res2.name.tolist().index("NH2")
            coor3 = np.copy(coor2)
            coor3[idx_nh1], coor3[idx_nh2] = coor2[idx_nh2], coor2[idx_nh1]
            diff = (coor1 - coor2).ravel()
            diff2 = (coor1 - coor3).ravel()
            return min(
                np.sqrt(3 * np.inner(diff, diff) / diff.size),
                np.sqrt(3 * np.inner(diff2, diff2) / diff2.size),
            )
        else:
            diff = (coor1 - coor2).ravel()
            return np.sqrt(3 * np.inner(diff, diff) / diff.size)

def get_chi1_chi2_rotamers(resname,rota_str):
    if rota_str=="OUTLIER":
        return rota_str
    if len(resname)==4:
        resname=resname[1:]
    if resname in ["MET","GLU","GLN","ARG","LYS"]:
        return rota_str[0:2]
    else:
        return rota_str
    

def get_rotamers(resol,run_num,ref_struct=False):
    if ref_struct:
        rotafile=open(base_dir+"7kr0_orig_rotalyze.out")
    else:
        rotafile=open(base_dir+"resol_"+resol+"/run_"+run_num+"_final/qfit_bic/rotalyze.out")
    qfit_rota_dict=dict()
    x=rotafile.readline()
    prev_resnum=-999
    first_res=True
    rota_array=[]
    for line in rotafile:
        if line.startswith("SUMMARY"):
            break
        lineparts=line[0:-1].split(":")
        resparts=lineparts[0].split()
        resnum=int(resparts[1])
        resname=resparts[2]
        if first_res:
            rotamers=get_chi1_chi2_rotamers(resname,lineparts[-1])
            rota_array.append(rotamers)
            first_res=False
        elif resnum!=prev_resnum:
            rotamer_set=sorted(set(rota_array))
            qfit_rota_dict[prev_resnum]=rotamer_set
            rota_array=[]
        rotamers=get_chi1_chi2_rotamers(resname,lineparts[-1])
        rota_array.append(rotamers)
        prev_resnum=resnum
    rotafile.close()
    rotamer_set=sorted(set(rota_array))
    qfit_rota_dict[prev_resnum]=rotamer_set
    return qfit_rota_dict

def get_altconf_rmsd(qfit_res,ref_res):
    rmsd_dict=dict()   
    # print(qfit_res.resn[0])    
    qfit_altlocs=sorted(list(set(qfit_res.altloc)))
    ref_altlocs = sorted(list(set(ref_res.altloc)))
    for qfit_altloc in qfit_altlocs:
        qfit_conf=qfit_res.extract("altloc",qfit_altloc)
        if qfit_conf.resn[0] in ['ALA','GLY']:
            qfit_conf_sidechain=qfit_conf.extract("name",['N','CA', 'C', 'O','CB'])
        else:
            qfit_conf_sidechain=qfit_conf.extract("name",['N','CA', 'C', 'O','H','HA','HA2','HA3'],"!=")
        if len(qfit_conf_sidechain.name)==0:
            continue
        rmsd_dict[qfit_altloc]=[]
        for ref_altloc in ref_altlocs:
            ref_conf=ref_res.extract("altloc",ref_altloc)
            if ref_conf.resn[0] in ['ALA','GLY']:
                ref_conf_sidechain=ref_conf.extract("name",['N','CA', 'C', 'O','CB'])
            else:
                ref_conf_sidechain=ref_conf.extract("name",['N','CA', 'C', 'O','H','HA','HA2','HA3'],"!=")
            if len(ref_conf_sidechain.name)==0:
                continue
            # rmsd=qfit_conf_sidechain.rmsd(ref_conf_sidechain)
            rmsd=my_rmsd(qfit_conf_sidechain,ref_conf_sidechain)
            rmsd_dict[qfit_altloc].append(rmsd)
    return rmsd_dict

def matching_rmsds(rmsd_dict):
    less_than_1_count_overall=0
    if len(rmsd_dict.keys())==0:
        return True
    if len(rmsd_dict.keys())>1:
        for altl in rmsd_dict:
            less_rmsd_found=False
            rmsd_array=rmsd_dict[altl]
            for i in rmsd_array:
                if i < 0.5:
                    less_rmsd_found=True
                    break
            if less_rmsd_found:
                less_than_1_count_overall+=1
        if less_than_1_count_overall >= len(rmsd_dict.keys()):
            return True
        else:
            return False
    else:
        rmsd_array=rmsd_dict['']
        if max(rmsd_array)<0.5:
            return True
        else:
            return False


def same_altconfs(ground_res,rotamer_dict):
    ground_resi=ground_res.resi[0]
    if ground_res.resn[0] in ['ALA','GLY']:
        number_of_rota=1
    else:
        number_of_rota=len(rotamer_dict[ground_resi])
        if number_of_rota > 1:
            return False
    rmsds=get_altconf_rmsd(ground_res,ground_res)
    for altconf in rmsds.keys():
        rmsd_array=rmsds[altconf]
        if max(rmsd_array)<=0.5:
            continue
        else:
            return False
    return True


def plot_correct_modeled_fraction():
    plt.figure()
    correct_df=pd.read_csv(base_dir+"correct_modeled_qfit_residues_fraction.csv")
    ax1=sns.lineplot(data=correct_df, x="resolution", y="overall_correct",ci='sd', linewidth=2)
    plt.xlim(0.8,3.0)
    ax1.set_xticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0]) 
    plt.ylim(0.0,1.0)
    plt.xlabel("Resolution ($\AA$)",fontsize=18)
    plt.ylabel("Fraction of correctly\nmodeled residues",fontsize=18)
    plt.savefig(base_dir+"Correct_modeled_residues_fraction.png",dpi=300,bbox_inches="tight")
    
def plot_multiconf_fraction():
    plt.figure()
    altconf_frac_df=pd.read_csv(base_dir+"multiconformer_residues_fraction.csv")
    ax1=sns.lineplot(data=altconf_frac_df, x="resolution", y="fraction of multiconformers",ci='sd', linewidth=2)
    plt.xlim(0.8,3.0)
    ax1.set_xticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0])
    plt.xlabel("Resolution ($\AA$)",fontsize=18)
    plt.ylabel("Fraction of residues that\nare multiconformer",fontsize=18)
    plt.savefig(base_dir+"fraction_of_multiconformer_residues.png",dpi=300,bbox_inches="tight")




r10s = range(8, 31)
number_of_runs=10
ref_rota_dict=get_rotamers("0","0",ref_struct=True)
opfile=open(base_dir+"correct_modeled_qfit_residues_fraction.csv","w")
opfile2=open(base_dir+"residue_level_all_details.csv","w")
opfile3=open(base_dir+"qfit_residue_level_classification.csv","w")
opfile4=open(base_dir+"overall_classification.csv","w")
opfile5=open(base_dir+"multiconformer_residues_fraction.csv","w")

opfile.write("resolution,run,overall_correct\n")
opfile2.write("resolution,run,residue,residue name,qfit altloc number,original altloc number,qfit_multiconf,ref_altconf,correct rmsd,correct rotamers,classification,rmsds,qfit_rotamers,ref_rotamers\n")
opfile3.write("resolution,run,residue,residue name,qfit altloc number,original altloc number,qfit_multiconf,ref_altconf,correct rmsd,correct rotamers,classification\n")
opfile4.write("resolution,run,tp,tn,fp,fn\n")
opfile5.write("resolution,run,fraction of multiconformers\n")

for r10 in r10s:
    resol = r10 * 0.1
    resol = str(round(resol,1))
    ref_struct=Structure.fromfile(base_dir+"xray_noise_modeling/p1_box_pdbs/7kr0_shift_bfac_"+resol+"_p1.pdb").reorder()
    for run in range(1,number_of_runs+1):
        fp=0
        tp=0
        fn=0
        tn=0
        multiconf_mismatch_count=0
        other_mismatch_count=0
        wrong_rmsd_count=0
        wrong_rotamer_count=0
        rmsd_correct_count=0
        rmsd_correct=False
        rota_correct=False
        overall_correct=False
        qfit_multiconf=False
        ref_multiconf=False
        qfit_multiconf_count=0
        qfit_rota_dict=get_rotamers(resol,str(run))
        qfit_struct=Structure.fromfile(base_dir+"resol_"+resol+"/run_"+str(run)+"_final/qfit_bic/7kr0_resol_"+resol+"_shake_noisy_run_"+str(run)+"_qFit.pdb").reorder()
        for ref_residue in (
                ref_struct.extract("record", "ATOM")  
                .extract("resn", "HOH", "!=") 
                .extract(
                    "name", "H", "!="
                )  
                .extract(
                    "e", "H", "!="
                )
        ).residue_groups:
            classification=""
            ref_altlocs = sorted(list(set(ref_residue.altloc)))
            resi = ref_residue.resi[0]
            orig_occ=len(ref_altlocs)
            qfit_residue=qfit_struct.extract("record","ATOM").extract("resi",resi).extract("name", "H", "!=").extract("e", "H", "!=")
            if len(qfit_residue.altloc)==0:
                opfile2.write(resol+","+str(run)+","+str(resi)+",,,,,,,,fn,,,\n")
                opfile3.write(resol+","+str(run)+","+str(resi)+",,,,,,,,fn\n")
                fn+=1
                continue
            qfit_altlocs=sorted(list(set(qfit_residue.altloc)))
            qfit_occ=len(qfit_altlocs)
            rmsds=get_altconf_rmsd(qfit_residue,ref_residue)

            qfit_occ_number=len(qfit_altlocs)
            if len(qfit_altlocs) > 1 and '' in qfit_altlocs:
                qfit_occ_number-=1
            orig_occ_number=len(ref_altlocs)
            if len(ref_altlocs) > 1 and '' in ref_altlocs:
                orig_occ_number-=1
            
            if qfit_occ_number>1:
                qfit_multiconf=not same_altconfs(qfit_residue,qfit_rota_dict)
            else:
                qfit_multiconf=False
            if orig_occ_number > 1:
                ref_multiconf=not same_altconfs(ref_residue,ref_rota_dict)
            else:
                ref_multiconf=False
            
            if qfit_multiconf: 
                qfit_multiconf_count+=1
                if ref_multiconf:
                    if matching_rmsds(rmsds):
                        rmsd_correct=True
                    else:
                        rmsd_correct=False
                        wrong_rmsd_count+=1
                    if qfit_residue.resn[0] in ["ALA","GLY"]:
                        rota_correct=True
                    else:
                        if qfit_rota_dict[resi] == ref_rota_dict[resi]:
                            rota_correct=True
                        else:
                            rota_correct=False
                            wrong_rotamer_count+=1
                    if rmsd_correct and rota_correct:
                        tp+=1
                        classification="tp"
                        overall_correct=True
                    else:
                        fp+=1
                        classification="fp"
                        overall_correct=False
                        other_mismatch_count+=1
                else:
                    fp+=1
                    classification="fp"
                    overall_correct=False
                    multiconf_mismatch_count+=1
            elif not qfit_multiconf:
                if not ref_multiconf:
                    if matching_rmsds(rmsds):
                        rmsd_correct=True
                    else:
                        rmsd_correct=False
                        wrong_rmsd_count+=1
                    if qfit_residue.resn[0] in ["ALA","GLY"]:
                        rota_correct=True
                    else:
                        if qfit_rota_dict[resi] == ref_rota_dict[resi]:
                            rota_correct=True
                        else:
                            rota_correct=False
                            wrong_rotamer_count+=1
                    if rmsd_correct and rota_correct:
                        tn+=1
                        classification="tn"
                        overall_correct=True
                    else:
                        fn+=1
                        classification="fn"
                        overall_correct=False
                        other_mismatch_count+=1
                else:
                    fn+=1
                    classification="fn"
                    overall_correct=False
                    multiconf_mismatch_count+=1
                
            
     
            opfile2.write(resol+","+str(run)+","+str(resi)+","+qfit_residue.resn[0]+","+str(qfit_occ_number)+","+str(orig_occ_number)+","+str(qfit_multiconf)+","+str(ref_multiconf)+","+str(rmsd_correct)+","+str(rota_correct)+","+classification+","+str(rmsds)+",")
            if resi not in qfit_rota_dict.keys():
                opfile2.write(",\n")
            else:
                opfile2.write(str(qfit_rota_dict[resi])+","+str(ref_rota_dict[resi])+"\n")    
            opfile3.write(resol+","+str(run)+","+str(resi)+","+qfit_residue.resn[0]+","+str(qfit_occ_number)+","+str(orig_occ_number)+","+str(qfit_multiconf)+","+str(ref_multiconf)+","+str(rmsd_correct)+","+str(rota_correct)+","+classification+"\n")

        correct_total=tp+tn
        correct_rate=correct_total/169.0
        incorrect_total=fp+fn
        incorrect_rate=incorrect_total/169.0

        if qfit_multiconf_count>0:
            tp_rate=tp/float(qfit_multiconf_count)
            fp_rate=fp/float(qfit_multiconf_count)
        else:
            tp_rate=0.00
            fp_rate=0.00
        qfit_single_conf_count=169-qfit_multiconf_count
        tn_rate=tn/float(qfit_single_conf_count)
        fn_rate=fn/float(qfit_single_conf_count)

        opfile.write(resol+","+str(run)+","+str(round(correct_rate,2))+"\n")
        
        opfile4.write(resol+","+str(run)+","+str(tp/169.0)+","+str(tn/169.0)+","+str(fp/169.0)+","+str(fn/169.0)+"\n")
        opfile5.write(resol+","+str(run)+","+str(qfit_multiconf_count/169.0)+"\n")
            
        
        
opfile.close()
opfile2.close()
opfile3.close()
opfile4.close()
opfile5.close()

plot_multiconf_fraction()
plot_correct_modeled_fraction()
