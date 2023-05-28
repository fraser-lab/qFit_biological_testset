#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 17:25:20 2023

@author: ashraya
"""
import argparse
import numpy as np
import os
from cctbx.maptbx.box import shift_and_box_model
from iotbx.data_manager import DataManager
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description='Generate noisy synthetic data at different resolutions given a PDB file')
    parser.add_argument('--base_dir', type=str, required=True, help='The base directory containing the PDB file. This is where output files and directories will be created')
    parser.add_argument('--pdb', type=str, required=True, help='name of pdb file, string before the .pdb extension')
    parser.add_argument('--model_resol','-r', required=True, type=float, help='resolution of the ground truth model')
    parser.add_argument('--resol_range_start','-rs', required=True, type=float, help='starting value of the resolution range to be evaluated')
    parser.add_argument('--resol_range_end','-re', required=True, type=float, help='ending value of the resolution range to be evaluated')
    parser.add_argument('--resol_step','-rstep', required=True, type=float, default=0.1, help='step size of the resolution range to be evaluated')
    parser.add_argument('--bfactor_scale', type=float, default=1.0, help='bfactor to be added for every 0.1A resolution decrease')
    parser.add_argument('--shift_bfac', action='store_true', help='shift bfactors, remove water, change symmetry to P1 and place in box')
    parser.add_argument('--shake_data', action='store_true', help='shake atom positions randomly using phenix.pdbtools and generate structure factors for shaken model')
    parser.add_argument('--shake_multiplier', type=float, default=0.2, help='multiplier * resolution will determine rmsd of shaking')
    parser.add_argument('--add_complex_noise', action='store_true', help='add noise in complex space')
    parser.add_argument('--noise_multiplier', type=float, default=0.5, help='Fnoise=Finput+(sqrt(Finput)*random gaussian*noise_multiplier*resolution)')
    parser.add_argument('--generate_single_conf', action='store_true', help='generate single conformer model from the b-factor shifted model')
    
    return parser.parse_args()

def shift_bfactors(resol_array, model_resol, bfactor_scale, pdb, base_dir):
    
    logfile = open(base_dir + pdb + "/noise_generator.log","a")
    logfile.write("Shifting B-factors...\n")
    
    if not os.path.exists(base_dir+pdb+"/"+pdb+".pdb"):
        print("Input pdb file does not exist")
        return
    file_name=base_dir+pdb+"/"+pdb+".pdb"
    for resol in resol_array:
        resol_str=str(round(resol,1))
        logfile.write(resol_str+" ")
        resol_diff = resol - model_resol
        bfact_to_add = (resol_diff * bfactor_scale) / 0.1
        out_filename=base_dir + pdb + "/" + pdb + "_shift_bfac_resol_"+resol_str+".pdb"
        cmd = "phenix.pdbtools %s output.file_name=%s convert_to_isotropic=True remove='water' shift_b_iso=%f" \
            % (file_name, out_filename, bfact_to_add)
        os.system(cmd)
    
    logfile.write("\nShifted bfactors\n")
    logfile.close()

def place_model_in_box(resol_array, pdb, base_dir):
    
    logfile = open(base_dir + pdb + "/noise_generator.log","a")
    logfile.write("Placing model in p1 box...\n")
    if not os.path.exists(base_dir + pdb + "/bfactor_shifted_pdbs/"):
            os.mkdir(base_dir + pdb + "/bfactor_shifted_pdbs")
    for resol in resol_array:
        resol_str=str(round(resol,1))
        logfile.write(resol_str+" ")
        dm = DataManager()
        dm.set_overwrite(True)
        file_name = base_dir + pdb + "/" + pdb + "_shift_bfac_resol_" + resol_str + ".pdb"
        model = dm.get_model(filename = file_name)
        model_p1 = shift_and_box_model(model)
        new_filename = base_dir + pdb + "/" + pdb + "_shift_bfac_resol_" + resol_str + "_p1.pdb"
        dm.write_model_file(model_p1, new_filename)
        shutil.move(base_dir + pdb + "/" + pdb + "_shift_bfac_resol_" + resol_str + "_p1.pdb", base_dir + pdb + "/bfactor_shifted_pdbs/" + pdb + "_shift_bfac_resol_" + resol_str + "_p1.pdb")
        os.remove(base_dir + pdb + "/" + pdb + "_shift_bfac_resol_" + resol_str + ".pdb")
    logfile.write("\nPlaced model in box\n")
    logfile.close()

def remove_altconfs(resol_array, pdb, base_dir):
    logfile = open(base_dir + pdb + "/noise_generator.log","a")
    logfile.write("Removing altconfs...\n")
    if not os.path.exists(base_dir + pdb + "/single_conformer_pdbs/"):
        os.mkdir(base_dir + pdb + "/single_conformer_pdbs") 
    
    for resol in resol_array:
        resol_str = str(round(resol,1))
        logfile.write(resol_str+" ")
        in_filename = base_dir + pdb + "/bfactor_shifted_pdbs/" + pdb + "_shift_bfac_resol_" + resol_str + "_p1.pdb"
        out_filename = base_dir + pdb + "/single_conformer_pdbs/" + pdb + "_shift_bfac_resol_" + resol_str + "_p1_one_conf.pdb"
        cmd = "phenix.pdbtools remove_alt_confs=True %s output.file_name=%s" \
            % (in_filename, out_filename)
        os.system(cmd)
        
    logfile.write("\nRemoved altconfs\n")
    logfile.close()
    
def shake_and_generate_data(resol_array, shake_multiplier, pdb, base_dir):
    
    logfile = open(base_dir + pdb + "/noise_generator.log","a")
    logfile.write("Shaking coordinates and generating data...\n")
    
    if not os.path.exists(base_dir + pdb + "/shaken_pdbs/"):
        os.mkdir(base_dir + pdb + "/shaken_pdbs") 
    if not os.path.exists(base_dir + pdb + "/shaken_maps/"):
        os.mkdir(base_dir + pdb + "/shaken_maps")
        
    for resol in resol_array:
        resol_str = str(round(resol,1))
        logfile.write(resol_str+" ")
        shake_scale = shake_multiplier * resol
        in_filename = base_dir + pdb + "/bfactor_shifted_pdbs/" + pdb + "_shift_bfac_resol_" + resol_str + "_p1.pdb"
        out_filename = base_dir + pdb + "/shaken_pdbs/" + pdb + "_resol_" + resol_str + "_shaken.pdb"
        shake_cmd = "phenix.pdbtools %s output.file_name=%s shake=%f" \
            % (in_filename, out_filename, shake_scale)
        os.system(shake_cmd)
        in_filename = base_dir + pdb + "/shaken_pdbs/" + pdb + "_resol_" + resol_str + "_shaken.pdb"
        out_filename = base_dir + pdb + "/shaken_maps/" + pdb + "_resol_" + resol_str + "_shaken.mtz"
        if os.path.exists(out_filename): #If the output file exists, phenix mistakenly takes it as reference reflection data and throws error
            os.remove(out_filename)
        data_cmd = "phenix.fmodel %s k_sol=0.4 b_sol=45 high_resolution=%f r_free_flags_fraction=0.05 output.file_name=%s" \
            % (in_filename, resol, out_filename)
        os.system(data_cmd)
    logfile.write("\nShaken and generated data\n")
    logfile.close()
        

def generate_complex_space_noise_script(resol_array, noise_multiplier, pdb, base_dir):
    if not os.path.exists(base_dir + pdb + "/complex_noise_maps/"):
        os.mkdir(base_dir + pdb + "/complex_noise_maps")
    opfile=open(base_dir + pdb + "/add_complex_space_noise.sh", "w")
    for resol in resol_array:
        resol_str = str(round(resol,1))
        noise_scale = noise_multiplier * resol
        noise_scale=str(round(noise_scale,1))
        opfile.write("echo \""+resol_str+"\"\n")
        opfile.write("sftools << eof\n")
        opfile.write("read shaken_maps/" + pdb + "_resol_" + resol_str + "_shaken.mtz\n")
        opfile.write("Y\n")
        opfile.write("calc col sqf = col FMODEL 0.5 **\n")
        opfile.write("calc col randf = col sqf ran_g " + noise_scale + " * *\n")
        opfile.write("calc F col FSIM = col FMODEL col randf +\n")
        opfile.write("calc Q col SIGFSIM = 0.1 col fmodel *\n")
        opfile.write("calc F col FWT = col FSIM\n") 
        opfile.write("calc P col PHWT = col PHIFMODEL\n")
        opfile.write("write complex_noise_maps/" + pdb + "_resol_" + resol_str + "_shake_noisy.mtz COL FSIM SIGFSIM FWT PHWT R-free-flags PHIFMODEL\n")
        opfile.write("quit\n")
        opfile.write("eof\n\n")
    opfile.close()

def add_sigf_column(resol_array, pdb, base_dir):
    opfile=open(base_dir + pdb + "/add_sigf_column.sh", "w")
    for resol in resol_array:
        resol_str = str(round(resol,1))
        opfile.write("echo \""+resol_str+"\"\n")
        opfile.write("sftools << eof\n")
        opfile.write("read shaken_maps/" + pdb + "_resol_" + resol_str + "_shaken.mtz\n")
        opfile.write("Y\n")
        opfile.write("calc F col FSIM = col FMODEL\n")
        opfile.write("calc Q col SIGFSIM = 0.1 col fmodel *\n")
        opfile.write("calc F col FWT = col FSIM\n")
        opfile.write("calc P col PHWT = col PHIFMODEL\n")
        opfile.write("write shaken_maps/" + pdb + "_resol_" + resol_str + "_shake_refine.mtz COL FSIM SIGFSIM FWT PHWT PHIFMODEL R-free-flags\n")
        opfile.write("quit\n")
        opfile.write("eof\n\n")
    opfile.close()

def main():

    args = parse_args()
    if not os.path.isdir(args.base_dir + args.pdb):
        print("base directory or pdb directory not found")
        return
    logfile = open(args.base_dir + args.pdb + "/noise_generator.log","w")
    resol_true_end = round(args.resol_range_end+args.resol_step,1)
    resol_array = np.arange(args.resol_range_start, resol_true_end, args.resol_step)
    logfile.write(str(resol_array) + "\n")
    logfile.close()
    if args.shift_bfac:
        shift_bfactors(resol_array, args.model_resol, args.bfactor_scale, args.pdb, args.base_dir)
        place_model_in_box(resol_array, args.pdb, args.base_dir)
    if args.generate_single_conf:
        remove_altconfs(resol_array, args.pdb, args.base_dir)
    if args.shake_data:
        shake_and_generate_data(resol_array, args.shake_multiplier, args.pdb, args.base_dir)
    if args.add_complex_noise:
        generate_complex_space_noise_script(resol_array, args.noise_multiplier, args.pdb, args.base_dir)
    logfile = open(args.base_dir + args.pdb + "/noise_generator.log","a")
    logfile.write("\nDone\n")
    logfile.close()
    
    
    


if __name__ == '__main__':
    main()