import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 5})
import pandas as pd
import numpy as np

def createHeatmap(data, title, xAxisLabels, yAxisLabels, outputFile):
        fig, ax = plt.subplots()
        im = ax.imshow(data)
        ax.figure.colorbar(im, ax=ax, orientation='horizontal')

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(xAxisLabels)))
        ax.set_yticks(np.arange(len(yAxisLabels)))

        # ... and label them with the respective list entries
        ax.set_xticklabels(xAxisLabels)
        ax.set_yticklabels(yAxisLabels)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(yAxisLabels)):
            for j in range(len(xAxisLabels)):
                text = ax.text(j, i, round(data[i, j], 2),
                            ha="center", va="center", color="w")

        ax.set_title(title)
        fig.tight_layout()
        plt.savefig(outputFile)



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

    yAxisLabels = range(firstResidueNum[group],firstResidueNum[group]+ligandLength[group])

    for pdbID in groups[group]:
        dir = "../"+group+"/"+pdbID+"/"

        singleMutationFile = pd.read_csv(dir+pdbID+"_rigidity_metric.csv")



        heatmapData = np.zeros([20, ligandLength[group]])

        heatmapDataCount = np.zeros([20, ligandLength[group]])

        for index, row in singleMutationFile.iterrows():
            mutation = row[1] # G2843T
            value = np.abs(row[4])

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

        heatmapData /= np.max(heatmapData)

        #This is the 2d array of values for the heatmap
        data = heatmapData.transpose()   



        createHeatmap(data, pdbID + " Rigidity Metric", aminoAcids, yAxisLabels, "output/single/"+pdbID+'_RA.png')

