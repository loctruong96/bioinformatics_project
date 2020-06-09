import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 5})
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

