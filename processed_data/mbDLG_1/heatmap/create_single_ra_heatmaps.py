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

bondDataFile = open("../ligandBonds.csv", "r")

bondData = {}

for line in bondDataFile:
    line = line.split(",")
    if(line[3].isnumeric()):
        if(line[0] not in bondData.keys()):
            bondData[line[0]] = {}
        if(line[1] not in bondData[line[0]].keys()):
            bondData[line[0]][line[1]] = {}
        if(line[2] not in bondData[line[0]][line[1]].keys()):
            bondData[line[0]][line[1]][line[2]] = [int(line[3]), int(line[4][:1])]

bondDataFile.close()

for group in groups.keys():


    for pdbID in groups[group]:

        yAxisLabels = []
        N = ligandLength[group]
        for i in range(N):
                branch = ligandBranch[group]
                residueNum = str(firstResidueNum[group] + i)
                ligandWT = f" ({ligand[group][i]})"
                yAxisLabels.append(branch+residueNum+ligandWT)

        for i in range(len(yAxisLabels)):
            bonds = bondData[pdbID]["WT"][yAxisLabels[i].split(" ")[0]]
            if(bonds[0]):
                yAxisLabels[i] = "["+yAxisLabels[i]+"]"
            if(bonds[1]):
                yAxisLabels[i] += "*"

        dir = "../"+group+"/"+pdbID+"/"

        singleMutationFile = pd.read_csv(dir+pdbID+"_rigidity_metric.csv")



        heatmapData = np.zeros([20, ligandLength[group]])

        heatmapDataCount = np.zeros([20, ligandLength[group]])

        axisData = np.zeros([20, ligandLength[group], 2])

        for index, row in singleMutationFile.iterrows():
            mutation = row[1] # G2843T
            value = row[2]

            aa = mutation[-1]

            #X coordinate
            aaNum = aminoAcids.index(aa)

            #Y coordinate
            resNum = int(mutation[1:-1])
            resNum -= firstResidueNum[group]
            
            axisData[aaNum, resNum] = bondData[pdbID][mutation][mutation[:-1]]
            heatmapData[aaNum, resNum] += value
            heatmapDataCount[aaNum, resNum] += 1
        

        
        heatmapDataCount[heatmapDataCount == 0] = 1

        heatmapData /= heatmapDataCount

        # heatmapData /= np.max(heatmapData)

        #This is the 2d array of values for the heatmap
        data = heatmapData.transpose()   



        createHeatmap(data, axisData.transpose(), pdbID + " Rigidity Metric", aminoAcids, yAxisLabels, "output/single/"+pdbID+'_RA.png')

