import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 5})
import numpy as np



def createHeatmap(data, axisData, title, xAxisLabels, yAxisLabels, outputFile):
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
                valString = round(data[i, j], 2)

                if(axisData[0,i,j]):
                        valString = "["+str(valString)+"]"
                if(axisData[1,i,j]):
                        valString = str(valString)+"*"

                if(data[i, j] == 0):
                        valString = "WT"
                        if(yAxisLabels[i][0] == "["):
                                valString = "["+str(valString)+"]"
                        if(yAxisLabels[i][-1] == "*"):
                                valString = str(valString)+"*"
                text = ax.text(j, i, valString,
                            ha="center", va="center", color="w")

        ax.set_title(title)
        fig.tight_layout()
        plt.savefig(outputFile, dpi=300)

