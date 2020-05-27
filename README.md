# bioinformatics_project

Call ./aggregateLogs.sh to regenerate all of the processed csv files found in ./processed_data/

or recreate them individually with

call ./aggregate_py/aggregate_raw_XXX_data.py PDBID
where XXX is either "em_step", "ra_cluster" or "ra_constraint" and PDBID is the id of the PDB that has output in ./preliminary_exp_output