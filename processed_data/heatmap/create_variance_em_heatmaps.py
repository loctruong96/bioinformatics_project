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

        ax.set_ylabel('Residue ID (WT)')
        ax.set_xlabel('Mutated Amino Acids')

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(yAxisLabels)):
            for j in range(len(xAxisLabels)):
                text = ax.text(j, i, int(data[i, j]),
                            ha="center", va="center", color="w")

        ax.set_title(title)
        fig.tight_layout()
        plt.savefig(outputFile)




groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"],
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

ligand = {
    "mbDLG_1": "SYLVTSV",
    "mbDLG_2": "YLVTSV",
    "mbDLG_3": "HETSV"
}

ligandBranch = {
    "mbDLG_1": "G",
    "mbDLG_2": "B",
    "mbDLG_3": "B"
}

ligandLength = {
    "mbDLG_1": 7,
    "mbDLG_2": 6,
    "mbDLG_3": 5
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

        doubleMutationFile = pd.read_csv(dir+pdbID+"_em_metric_2.csv")



        heatmapDataAll = np.zeros([20, ligandLength[group], 19*ligandLength[group]])
        
        heatmapDataCurIndex = np.zeros([20, ligandLength[group]]).astype(np.int32)


        for index, row in doubleMutationFile.iterrows():
            mutation = row[1] # G2842A.G2843T
            mutations = mutation.split(".")
            
            value = np.abs(row[3])

            for mut in mutations:
                #mut = G2842A
                aa = mut[-1]

                #X coordinate
                aaNum = aminoAcids.index(aa)

                #Y coordinate
                resNum = int(mut[1:-1])
                resNum -= firstResidueNum[group]
                
                heatmapDataAll[aaNum, resNum, int(heatmapDataCurIndex[aaNum, resNum])] += value
                heatmapDataCurIndex[aaNum, resNum] += 1

        heatmapVariances = np.var(heatmapDataAll, axis=2)
        heatmapData = np.zeros([20, ligandLength[group]])

        for i in range(20):
            for j in range(ligandLength[group]):
                std = np.std(heatmapDataAll[i,j])
                mean = np.mean(heatmapDataAll[i,j])
                for k in range(19*ligandLength[group]):
                    if(heatmapDataAll[i,j,k] - mean > 1*std):
                        heatmapData[i,j] += 1


        #This is the 2d array of values for the heatmap
        data = heatmapData.transpose()   



        createHeatmap(data, pdbID + " EM Outliers by 1 Standard Deviation", aminoAcids, yAxisLabels, "output/variance/"+pdbID+'_EM.png')
