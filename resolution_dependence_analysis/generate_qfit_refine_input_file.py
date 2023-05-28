#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 01:27:24 2023

@author: ashraya
"""


import numpy as np
import argparse
def parse_args():
    parser = argparse.ArgumentParser(description='Create a file that lists pdbs to refine and the the synthetic resolutions')
    parser.add_argument('--pdb', type=str, required=True, help='Tname of pdb file, string before the .pdb extension')
    parser.add_argument('--base_dir', type=str, required=True, help='The base directory containing the PDB file. This is where output files and directories will be created')
    parser.add_argument('--resol_range_start','-rs', required=True, type=float, help='starting value of the resolution range to be evaluated')
    parser.add_argument('--resol_range_end','-re', required=True, type=float, help='ending value of the resolution range to be evaluated')
    parser.add_argument('--resol_step','-rstep', required=True, type=float, default=0.1, help='step size of the resolution range to be evaluated')
    return parser.parse_args()

def main():
    args = parse_args()
    opfile=open(args.base_dir + args.pdb + "/qfit_input_list.txt","w")
    resol_true_end = round(args.resol_range_end+args.resol_step,1)
    resol_array = np.arange(args.resol_range_start, resol_true_end, args.resol_step)
    for resol in resol_array:
        opfile.write(args.pdb + "," + str(round(resol,1)) + "\n")
    opfile.close()

if __name__ == '__main__':
    main()