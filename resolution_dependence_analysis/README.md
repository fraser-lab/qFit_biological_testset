# qFit_biological_testset

Prerequisites: PHENIX, CCP4 and qFit
Instructions to run the qFit resolution dependence analysis:

1. Have the PDB file (7kr0.pdb for the paper) of the ground truth structure of interest in a base directory (say base_dir)
2. Run make_some_noise.py with appropriate command line arguments. For the paper, following command was run:
`cctbx.python make_some_noise.py --base_dir=<base_dir> --pdb=7kr0 -r=0.77 -rs=0.8 -re=3.0 -rstep=0.1`
	This will create a folder with the name of pdb file and some files and folders within the pdb directory
3. Run the *add_sigf_column.sh* script generated inside the pdb folder. This adds the SIGFSIM column to the mtz file which is required for refinement
4. Run the *add_complex_space_noise.sh* script generated inside the pdb folder. This adds complex noise to structure factors of the shaken mtz files
5. Run the *generate_qfit_refine_input_file.py* script as follows:
`python generate_qfit_refine_input_file.py --pdb=7kr0 --base_dir=<base_dir> -rs=0.8 -re=3.0 -rstep=0.1`
	This creates a file called *qfit_input_list.txt*
6. Run the *pre_qfit_refine_parallel.sge* script. This is script (and other .sge scripts) has been written to work with a multi-cpu workstation with SGE.
7. Run the *qfit_qsub_run.sge* script
8. Run the *post_qfit_qsub_run.sge* script
	After all these steps, there will be separate directories created for each resolution and within that the final qFit model will be named *7kr0_resol_RESOL_shake_noisy_qFit.pdb*, where RESOL is the resolution of the dataset.
9. Repeat steps 2-8 for multiple runs (10 runs in the paper). Make sure to change base directory for each run to prevent overwriting

Scripts for reproducing figures:
Before running scripts for generating figures, the rotamer states of residues need to be calculated (using *run_rotalyze.sh*)
1. *evaluate_correct_modeled_conformers.py* - This reproduces figure 4A and 4B
2. *mapq_calculate_qscore.sh* - This script calculates Q-score for every atom in the input structure using mapq package along with Chimera. This done for the ground truth models and qFit output models
3. *compare_qscore.py* - This reproduces figure 4D

