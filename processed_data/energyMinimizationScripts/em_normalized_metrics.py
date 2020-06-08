import pandas as pd

groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"]
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

for group in groups.keys():
    for pdb in group:
        dir = group+"/"+pdbID+"/
        with open(+pdbID+"_normalized_em_metrics.csv") as output:
