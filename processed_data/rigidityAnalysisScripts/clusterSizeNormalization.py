import sys
import csv
import os


def process(pdbID):
    dataFile = "../"+pdbID+"/"+pdbID+"_raw_ra_cluster_data.csv"
    baselineFile = "../"+pdbID+"/base_"+pdbID+"_raw_ra_cluster_data.csv"


    baselineData = {}
    with open(baselineFile) as csvfile:
        baselineList = list(csv.reader(csvfile))
        for i in range(1, len(baselineList)):
            baselineData[baselineList[i][2]] = baselineList[i][3]

    out = open("../"+pdbID+"/"+pdbID+"_normalized_clusters.csv","w+")
    writer = csv.writer(out, delimiter=",")


    with open(dataFile, "r") as inFile:
        reader = csv.reader(inFile)
        header = next(reader)
        writer.writerow(header)
        for row in reader:
            if row[2] in baselineData.keys(): 
                #print("baseline for {0}: {1}".format(row[2], baselineData[row[2]]))
                row[3] = int(row[3]) - int(baselineData[row[2]])
            if row[3] > 0:
                writer.writerow(row)

    out.close()

    print(baselineData)




if __name__ == '__main__':
    assert len(sys.argv) > 1
    process(sys.argv[1].upper())



