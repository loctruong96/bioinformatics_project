This folder contain energy minimization and rigidity analysis baselines for '2DM8', 'MBDLG', and '3RL7' pdb files. 

Note that the mutation being performed on each of these pdb files was G2837S. Meaning, residue 2837 on chain G is being mutated from S -> S ... aka ... no mutation was done.

The redundancy check in promute source code (which prevent S -> S) was disabled to generate these baselines (in case we want to regenerate the baselines for comparison). 

Similar to preliminary_exp_output's structure, multiMutant_out contains output from multiMutant pipeline and additional_EM_logs folder contains EM logs that are not included in multiMutant pipeline.


