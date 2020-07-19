import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm, matplotlib.colors
plt.rcParams.update({'font.size': 5})
import numpy as np

import copy
from createHeatmap import createHeatmap

'''
def createHeatmap2(data, min, title, xAxisLabels, yAxisLabels, outputFile):
        fig, ax = plt.subplots()

        my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default
        my_cmap.set_bad((0,0,0))

        # im = ax.imshow(data,norm=matplotlib.colors.SymLogNorm(linthresh=min,vmin=data.min(),vmax=data.max()), cmap=my_cmap)
        im = ax.imshow(data)
        ax.figure.colorbar(im, ax=ax, orientation='horizontal')

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(xAxisLabels)))
        ax.set_yticks(np.arange(len(yAxisLabels)))

        # ... and label them with the respective list entries
        ax.set_xticklabels(xAxisLabels)
        ax.set_yticklabels(yAxisLabels)

        ax.set_ylabel('Residue ID (WT)')
        ax.set_xlabel('Mutated Amino Acids')

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(yAxisLabels)):
            for j in range(len(xAxisLabels)):
                valString = "0"
                if data[i, j] > 0:
                    valString = "{:.0e}".format(data[i, j])
                text = ax.text(j, i, valString,
                            ha="center", va="center", color="w")

        ax.set_title(title)
        fig.tight_layout()
        plt.savefig(outputFile, dpi=300)
'''


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

        emMetric = pd.read_csv(dir+pdbID+"_em_metric.csv")
        raMetric = pd.read_csv(dir+pdbID+"_rigidity_metric.csv")      

        emheatMap = np.zeros([20, ligandLength[group]])
        raheatMap = np.zeros([20, ligandLength[group]])
        axisData = np.zeros([20, ligandLength[group], 2])

        for index, row in emMetric.iterrows():
            mutation = row[1] # G2843T
            value = row[3]

            aa = mutation[-1]

            #X coordinate
            aaNum = aminoAcids.index(aa)

            #Y coordinate
            resNum = int(mutation[1:-1])
            resNum -= firstResidueNum[group]

            axisData[aaNum, resNum] = bondData[pdbID][mutation][mutation[:-1]]
            emheatMap[aaNum, resNum] += value


        for index, row in raMetric.iterrows():
            mutation = row[1] # G2843T
            value = np.abs(row[4])

            aa = mutation[-1]

            #X coordinate
            aaNum = aminoAcids.index(aa)

            #Y coordinate
            resNum = int(mutation[1:-1])
            resNum -= firstResidueNum[group]
            
            raheatMap[aaNum, resNum] += value
        
        raheatMap /= raheatMap.max()


        combinedHeatMap = emheatMap * raheatMap

        data = combinedHeatMap.transpose()


        min = data[data > 0].min()

        createHeatmap(data, axisData.transpose(), pdbID+" EM Metric * RA Metric", aminoAcids, yAxisLabels, "output/combined/"+pdbID+'_combined.png')

