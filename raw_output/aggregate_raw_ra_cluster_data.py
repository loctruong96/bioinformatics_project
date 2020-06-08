import sys
import os
import csv


baseline = ""
baseFileName = ""
if(len(sys.argv) > 1):
    if sys.argv[1] == "-b":
        baseline += "baselines/"
        baseFileName = "base_"

groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"],
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

for group in groups.keys():
    for pdbID in groups[group]:
        dir = group+"_out/"+baseline+"multiMutant_out/"




        outFilename = baseFileName+pdbID+"_raw_ra_cluster_data.csv"

        tmp = os.listdir(dir)
        mutName = ""
        for folder in tmp:
            if folder[:len(pdbID)] == pdbID:
                mutName = folder
                break
        dir+=mutName
        #example dir: "/preliminary_exp_output/multiMutant_out/2DM8G2837:2843_out/"

        #Output configuration

        columns = ["PDBID", "MUTATION", "CLUSTER_SIZE", "COUNT"]

        outDir = "../processed_data/"+group+"/"+pdbID

        if not os.path.exists("../processed_data"):
            os.makedirs("../processed_data")
        if not os.path.exists("../processed_data/"+group):
            os.makedirs("../processed_data/"+group)
        if not os.path.exists(outDir):
            os.makedirs(outDir)


        out = open(outDir+"/"+outFilename,"w+")
        writer = csv.writer(out, delimiter=",")
        writer.writerow(columns)

        numRows = 0

        mutations = os.listdir(dir)

        foundFiles = 0
        for mutation in mutations:
            fileNameRoot = "{0}/{1}/{2}.all.ra.out_RA/{2}.all.processed.pdb/user/{2}.all.processed.pdb".format(dir,mutation,pdbID)
            if(len(pdbID) > 4):
                fileNameRoot = "{0}/{1}/{3}.all.ra.out_RA/{3}.all.processed.pdb/user/{3}.all.processed.pdb".format(dir,mutation,pdbID,pdbID[:4])
            #filename = fileNameRoot+"_MetricsPDB.txt"
            filename = fileNameRoot+"_postPG_MetricsBBH.txt"
            if not os.path.exists(filename):
                fileNameRoot = "{0}/{1}/{2}.all.ra.out_RA/{2}.all.ra.out/{2}.all.processed.pdb/user/{2}.all.processed.pdb".format(dir,mutation,pdbID)
                if(len(pdbID) > 4):
                    fileNameRoot = "{0}/{1}/{3}.all.ra.out_RA/{3}.all.ra.out/{3}.all.processed.pdb/user/{3}.all.processed.pdb".format(dir,mutation,pdbID,pdbID[:4])
            filename = fileNameRoot+"_postPG_MetricsBBH.txt"
            if os.path.exists(filename):
                log = open(filename, "r") 
                splitMut = mutation.split(".")
                
                # mutationStr = splitMut[1]+"."+splitMut[2]
                mutationStr = splitMut[1]
                
                foundFiles += 1

                for line in log:
                    if(line[:4] == "size"):
                        split = line[5:-1].split(": ")
                        if(split[0].isdigit()):
                            data = [pdbID, mutationStr, split[0], split[1]]
                            writer.writerow(data)
                            numRows += 1
                log.close()

        out.close()
        print "Aggregated {0} logs into {1}, {2} rows".format(foundFiles, outFilename, numRows)