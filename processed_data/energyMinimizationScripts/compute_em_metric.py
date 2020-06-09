import pandas as pd
import numpy as np
from scipy.special import softmax
from scipy import stats
import numpy as np

groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"],
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

for group in groups.keys():
    for pdbID in groups[group]:
        for i in range(2):
            dir = "../"+group+"/"+pdbID+"/"
            metricsFile = dir+pdbID+"_refined_em_metrics"
            outputFile = dir+pdbID+"_em_metric"
            if(i == 1):
                metricsFile += "_2"
                outputFile += "_2"
            
            allMetrics = pd.read_csv(metricsFile+".csv")
            

            D_WT10_95 = allMetrics.D_WT10_95
           
            never95 = allMetrics[allMetrics.WT10_95 < 0].index.tolist()
            WT10_95 = allMetrics.WT10_95
            tmp = D_WT10_95

            for index in never95:
                D_WT10_95[index] = 0


            sm = softmax(D_WT10_95.abs())
            std = np.std(D_WT10_95.abs())

            #norm = D_WT10_95.abs() / std
            norm = D_WT10_95.abs()
            norm = norm / np.max(norm)

            # print("std")
            # print(std)

            for index in never95:
                sm[index] = 1
                norm[index] = 1

            # z = stats.zscore(sm)

            d = {
                "PDBID": allMetrics.PDBID,
                "MUTATION": allMetrics.MUTATION,
                "EM_METRIC": D_WT10_95,
                "NORMALIZED_EM_METRIC": norm,
                "SOFTMAX_EM_METRIC": sm
                # "Z_SCORE": z,
            }

            df = pd.DataFrame(data=d)

            df.to_csv(outputFile+".csv", index=False, 
            columns= ["PDBID", 
                "MUTATION", 
                "EM_METRIC", 
                "NORMALIZED_EM_METRIC", 
                "SOFTMAX_EM_METRIC",
                # "Z_SCORE"
            ])
            