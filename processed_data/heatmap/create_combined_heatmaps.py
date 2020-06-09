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
