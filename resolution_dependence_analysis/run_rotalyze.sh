#!/bin/bash

#This script runs phenix.rotalysze to calculate the rotamer angles and states
base_dir='/Users/ashraya/qfit_summit/summit_2/macdomain_maps/'

for resol in 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0
do
	echo $resol
	for run in 1 2 3 4 5 6 7 8 9 10
	do
		phenix.rotalyze ${base_dir}resol_${resol}/run_${run}_final/qfit_bic/7kr0_resol_${resol}_shake_noisy_run_${run}_qFit.pdb >> ${base_dir}resol_${resol}/run_${run}_final/qfit_bic/rotalyze.out
	done
done