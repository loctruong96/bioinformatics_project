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

        singleMutationFile = pd.read_csv(dir+pdbID+"_em_metric.csv")


        heatmapData = np.zeros([20, ligandLength[group]])

        heatmapDataCount = np.zeros([20, ligandLength[group]])

        for index, row in singleMutationFile.iterrows():
            mutation = row[1] # G2843T
            value = row[3]

            aa = mutation[-1]

            #X coordinate
            aaNum = aminoAcids.index(aa)

            #Y coordinate
            resNum = int(mutation[1:-1])
            resNum -= firstResidueNum[group]

            heatmapData[aaNum, resNum] += value
            heatmapDataCount[aaNum, resNum] += 1

        heatmapDataCount[heatmapDataCount == 0] = 1

        heatmapData /= heatmapDataCount

        #This is the 2d array of values for the heatmap
        data = heatmapData.transpose()

        createHeatmap(data, pdbID+" Energy Minimization Metric", aminoAcids, yAxisLabels, "output/single/"+pdbID+'_EM.png')
