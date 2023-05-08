#________________________________________________SET PATHS________________________________________________#
source /wynton/group/fraser/swankowicz/phenix-installer-1.20.1-4487-intel-linux-2.6-x86_64-centos6/phenix-1.20.1-4487/phenix_env.sh
export PATH="/wynton/home/fraserlab/swankowicz/anaconda3/bin:$PATH"
source activate qfit
#which python
export PHENIX_OVERWRITE_ALL=true

#________________________________________________RUN QFIT________________________________________________#
 PDB_dir=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)
 echo ${PDB_dir}
 cd ${base_dir}/${PDB_dir}
 pwd
 PDB=$(echo ${PDB_dir:0:4})
 echo ${PDB}
 mkdir ${category}
 cd ${category}
 pwd
 #cp ../${PDB}.mtz .
 #cp ../composite_omit_map.mtz .
 #cp ../${PDB}.pdb .
 #qfit_protein composite_omit_map.mtz -l 2FOFCWT,PH2FOFCWT ${PDB}.pdb -p 8
 qfit_final_refine_xray.sh ${PDB}.mtz multiconformer_model2.pdb
