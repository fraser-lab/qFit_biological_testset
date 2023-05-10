#________________PDB INFO__________________________________
PDB_file=/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/list_of_pdbs_final.txt
PDB_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/'
output_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/output_data/nconfs/'

category='nconf'

PDB=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)



#______________________PDB STATS_______________________

   qfit_RMSF.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${PDB}
   calc_occ.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb ${PDB}
   b_fac=$(b_factor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb=${PDB})
   phenix.rotalyze model=${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb outliers_only=False > ${output_dir}/${PDB}_rotamer_output.txt
   python /wynton/group/fraser/swankowicz/script/refine_log_parse.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log ${PDB_dir}/${PDB}/${PDB}_single_001.log ${PDB}
   mmtbx.validation_summary ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_validation.txt
