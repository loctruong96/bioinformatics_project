import sys
import csv
import json
import os
import clusterSizeNormalization


def process(pdbID):

    dataFile = "../"+pdbID+"/"+pdbID+"_normalized_clusters.csv"

    if not os.path.isfile(dataFile):
        clusterSizeNormalization.process(pdbID)

    outliers = {}

    with open(dataFile, "r") as inFile:
        reader = csv.reader(inFile)
        header = next(reader)
        for row in reader:
            if row[1] not in outliers.keys(): 
                outliers[row[1]] = {}
                outliers[row[1]]["count"] = 0
                outliers[row[1]]["clusters"] = []
            if int(row[2]) > 100:
                outliers[row[1]]["count"] += 1
                outliers[row[1]]["clusters"].append([row[2],row[3]])

    with open("../"+pdbID+"/"+pdbID+"_cluster_outliers.json", 'w') as json_file:
        json.dump(outliers, json_file)

    out = open("../"+pdbID+"/"+pdbID+"_cluster_outlier_count.csv","w+")
    writer = csv.writer(out, delimiter=",")
    writer.writerow(["MUTATION","OUTLIERS"])

    for mutant in outliers:
        count = outliers[mutant]["count"]
        writer.writerow([mutant,count])

    out.close()

    for mutant in outliers:
        count = outliers[mutant]["count"]
        if count > 1:
            print(mutant+": " + str(count))


if __name__ == '__main__':
    assert len(sys.argv) > 1
    process(sys.argv[1].upper())