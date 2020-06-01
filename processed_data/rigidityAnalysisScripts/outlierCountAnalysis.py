import sys
import csv
import os


def process(pdbID): 
    dataFile = "../"+pdbID+"/"+pdbID+"_cluster_outlier_count.csv"

    if not os.path.isfile(dataFile):
        clusterSizeNormalization.process(pdbID)


    with open(dataFile) as csvfile:
        outliers = list(csv.reader(csvfile))[1:]


    #By Residue
    residueCounts = {}
    for count in outliers:
        residue = count[0][1:-1]
        if residue not in residueCounts.keys():
            residueCounts[residue] = 0
        residueCounts[residue] += int(count[1])

    #By AA
    aaCounts = {}
    for count in outliers:
        aa = count[0][-1]
        if aa not in aaCounts.keys():
            aaCounts[aa] = 0
        aaCounts[aa] += int(count[1])


    #By AA group

    posCharged = ["R", "H", "K"]
    negCharged = ["D", "E"]
    polarUncharged = ["S", "T", "N", "Q"]
    special = ["C", "G", "P"]
    hydrophobic = ["A", "V", "I", "L", "M", "F", "Y", "W"]

    aaGroups = [posCharged, negCharged, polarUncharged, special, hydrophobic]
    groupCounts = [0.0,0.0,0.0,0.0,0.0]
    for count in outliers:
        aa = count[0][-1]
        for i in range(len(aaGroups)):
            if aaGroups[i].count(aa) > 0:
                groupCounts[i] += int(count[1])
                break
    for i in range(len(aaGroups)):
        groupCounts[i] = groupCounts[i] / len(aaGroups[i])
    # print(residueCounts)
    # print(aaCounts)
    # print(groupCounts)


    out = open("../"+pdbID+"/"+pdbID+"_cluster_outlier_analysis.txt","w+")
    out.write(pdbID + "Cluster Outlier Analysis\n")
    out.write("\n\nBy Residue:\n")
    out.write("===============\n")
    while len(residueCounts.keys()) > 0:
        max = -1
        res = ""
        for residue in residueCounts.keys():
            if residueCounts[residue] > max:
                res = residue
                max = residueCounts[residue]
        residueCounts.pop(res)
        out.write(res + ": " + str(max) + "\n")



    # for residue in residueCounts.keys():
    #     out.write(residue + ": " + str(residueCounts[residue]) + "\n")
    out.write("\n\nBy Amino Acid\n")
    out.write("===============\n")
    while len(aaCounts.keys()) > 0:
        max = -1
        a = ""
        for aa in aaCounts.keys():
            if aaCounts[aa] > max:
                a = aa
                max = aaCounts[aa]
        aaCounts.pop(a, None)
        out.write(a + ": " + str(max) + "\n")


    # for aa in aaCounts.keys():
    #     out.write(aa + ": " + str(aaCounts[aa]) + "\n")

    out.write("\n\nBy Amino Acid Group (Count / Size of Group)\n")
    out.write("===============\n")
    out.write("Positively Charged: " + str(groupCounts[0])+"\n")
    out.write("Negatively Charged: " + str(groupCounts[1])+"\n")
    out.write("Polar Uncharged   : " + str(groupCounts[2])+"\n")
    out.write("Special Cases     : " + str(groupCounts[3])+"\n")
    out.write("Hydrophobic       : " + str(groupCounts[4])+"\n")

    out.close()

if __name__ == '__main__':
    assert len(sys.argv) > 1
    process(sys.argv[1].upper())



