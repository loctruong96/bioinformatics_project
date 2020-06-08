import pandas as pd
import numpy as np
from scipy import stats

groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"],
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

for group in groups.keys():
    for pdbID in groups[group]:
        for i in range(2):
            dir = "../"+group+"/"+pdbID+"/"
            baselineFile = dir+"base_"+pdbID+"_em_metrics"
            metricsFile = dir+pdbID+"_em_metrics"
            outputFile = dir+pdbID+"_refined_em_metrics"
            if(i == 1):
                metricsFile += "_2"
                outputFile += "_2"
            
            baseline = pd.read_csv(baselineFile+".csv")

            allMetrics = pd.read_csv(metricsFile+".csv")
            

            mt10_95 = allMetrics.STEPS_MT10_95
            wt10_95 = allMetrics.STEPS_WT10_95

            norm_mt10_95 = mt10_95 - baseline.STEPS_MT10_95[0]
            norm_wt10_95 = wt10_95 - baseline.STEPS_WT10_95[0]

            mtZ = stats.zscore(mt10_95)
            wtZ = stats.zscore(wt10_95)

            nmtZ = stats.zscore(norm_mt10_95)
            nwtZ = stats.zscore(norm_wt10_95)

            d = {
                "PDBID": allMetrics.PDBID,
                "MUTATION": allMetrics.MUTATION,
                "MT10_95": mt10_95,
                "Z_MT10_95": mtZ,
                "D_MT10_95": norm_mt10_95,
                "Z_D_MT10_95": nmtZ,
                "WT10_95": wt10_95,
                "Z_WT10_95": wtZ,
                "D_WT10_95": norm_wt10_95,
                "Z_D_WT10_95": nwtZ,
                "BELOW_WT": allMetrics.STEPS_WT_100
            }

            df = pd.DataFrame(data=d)

            df.to_csv(outputFile+".csv", index=False, 
            columns= ["PDBID", "MUTATION", "MT10_95", "Z_MT10_95", "D_MT10_95", "Z_D_MT10_95", "WT10_95", "Z_WT10_95", "D_WT10_95", "Z_D_WT10_95", "BELOW_WT"])
            