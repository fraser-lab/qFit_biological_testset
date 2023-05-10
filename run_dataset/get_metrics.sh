#!/bin/bash
#$ -l h_vmem=4G
#$ -l mem_free=4G
#$ -t 1-174
#$ -l h_rt=2:00:00
'''
This will run all of the metrics we used for the qFit paper on an SGE server. 
qFit_RMSF.py and b_factor.py are loaded with qFit (https://github.com/ExcitedStates/qfit-3.0)
phenix.rotalyze and mmtbx.validation_summary are loaded with the Phenix program (https://phenix-online.org/)
refine_log_parse.py is located in this repository.
'''

#__________________SOURCE PHENIX/QFIT________________________________________________#
source phenix_env.sh
export PATH="/wynton/home/fraserlab/swankowicz/anaconda3/bin:$PATH"
source activate qfit
export PHENIX_OVERWRITE_ALL=true

#________________PDB INFO__________________________________
PDB_file=./expanded_dataset/list_of_pdbs_final.txt
PDB_dir='./expanded_dataset/'
output_dir='./output_data/nconfs/'

category='nconf'

PDB=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)


#______________________PDB STATS_______________________

qfit_RMSF.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${PDB}
b_fac=$(b_factor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb=${PDB})
phenix.rotalyze model=${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb outliers_only=False > ${output_dir}/${PDB}_rotamer_output.txt
mmtbx.validation_summary ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_validation.txt
python refine_log_parse.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log ${PDB_dir}/${PDB}/${PDB}_single_001.log ${PDB}

