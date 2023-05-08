#!/bin/bash
#$ -l h_vmem=4G
#$ -l mem_free=4G
#$ -t 1-144
#$ -l h_rt=24:00:00
#$ -R yes
#$ -pe smp 8

#________________________________________________INPUTS________________________________________________#
PDB_file=/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/list_of_pdbs_final.txt
base_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset'
category='test'
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


#________________________________________________SET PATHS________________________________________________#
#source phenix
source phenix_env.sh

#activate qfit
source activate qfit
export PHENIX_OVERWRITE_ALL=true

#________________________________________________RUN QFIT________________________________________________#
 PDB=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)
 echo ${PDB}
 cd ${base_dir}/${PDB_dir}
 
 mkdir ${category}
 cd ${category}
 
 cp ../${PDB}.mtz .
 cp ../composite_omit_map.mtz .
 cp ../${PDB}.updated_refine_001.pdb .
 
 qfit_protein composite_omit_map.mtz -l 2FOFCWT,PH2FOFCWT ${PDB}.updated_refine_001.pdb -p 8
 qfit_final_refine_xray.sh ${PDB}.mtz multiconformer_model2.pdb
