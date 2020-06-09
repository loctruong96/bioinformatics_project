from createHeatmap import createHeatmap
import pandas as pd
import numpy as np

groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"],
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

ligandLength = {
    "mbDLG_1": 7,
    "mbDLG_2": 6,
    "mbDLG_3": 5
}

ligandBranch = {
    "mbDLG_1": "G",
    "mbDLG_2": "B",
    "mbDLG_3": "B"
}

ligand = {
    "mbDLG_1": "SYLVTSV",
    "mbDLG_2": "YLVTSV",
    "mbDLG_3": "HETSV"
}

firstResidueNum = {
    "mbDLG_1": 2837,
    "mbDLG_2": 2838,
    "mbDLG_3": 513
}

aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
 				    "M", "F", "P", "S", "T", "W", "Y", "V"]

for group in groups.keys():

    yAxisLabels = []
    N = ligandLength[group]
    for i in range(N):
            branch = ligandBranch[group]
            residueNum = str(firstResidueNum[group] + i)
            ligandWT = f" ({ligand[group][i]})"
            yAxisLabels.append(branch+residueNum+ligandWT)

    for pdbID in groups[group]:
        dir = "../"+group+"/"+pdbID+"/"

        emMetric = pd.read_csv(dir+pdbID+"_em_metric.csv")
        raMetric = pd.read_csv(dir+pdbID+"_rigidity_metric.csv")
        doubleMutationFile = pd.read_csv(dir+pdbID+"_em_metric_2.csv")

        emheatMap = np.zeros([20, ligandLength[group]])
        raheatMap = np.zeros([20, ligandLength[group]])
        for index, row in emMetric.iterrows():
            mutation = row[1] # G2843T
            value = row[3]

            aa = mutation[-1]

            #X coordinate
            aaNum = aminoAcids.index(aa)

            #Y coordinate
            resNum = int(mutation[1:-1])
            resNum -= firstResidueNum[group]
            
            emheatMap[aaNum, resNum] += value

        for index, row in raMetric.iterrows():
            mutation = row[1] # G2843T
            value = row[3]

            aa = mutation[-1]

            #X coordinate
            aaNum = aminoAcids.index(aa)

            #Y coordinate
            resNum = int(mutation[1:-1])
            resNum -= firstResidueNum[group]
            
            raheatMap[aaNum, resNum] += value
        
        raheatMap /= np.max(raheatMap)

        #This is the 2d array of values for the heatmap

        combinedHeatMap = emheatMap * raheatMap

        print(pdbID)
        print(np.max(raheatMap))
        print(np.max(emheatMap))
        print(np.max(combinedHeatMap))
        print("\n")
        data = combinedHeatMap.transpose()   

        createHeatmap(data, pdbID+" Energy Minimization Metric", aminoAcids, yAxisLabels, "output/combined/"+pdbID+'_combined.png')
