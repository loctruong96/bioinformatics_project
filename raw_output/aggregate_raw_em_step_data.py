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
        dir = group+"_out/"+baseline+"additional_EM_logs/"+pdbID+"_EM_logs/"

        outFileName = baseFileName+pdbID+"_raw_em_step_data.csv"

        #Output configuration

        columns = ["PDBID", "MUTATION", "TS", "POTENTIAL"]

        outDir = "../processed_data/"+group+"/"+pdbID

        if not os.path.exists("../processed_data"):
            os.makedirs("../processed_data")
        if not os.path.exists("../processed_data/"+group):
            os.makedirs("../processed_data/"+group)
        if not os.path.exists(outDir):
            os.makedirs(outDir)


        out = open(outDir+"/"+outFileName,"w+")
        writer = csv.writer(out, delimiter=",")
        writer.writerow(columns)

        numRows = 0

        mutations = os.listdir(dir)

        foundFiles = 0

        for mutation in mutations:
            filename = dir + mutation + "/" + mutation + "_min.log"
            if os.path.exists(filename):
                log = open(filename, "r") 
                splitMut = mutation.split(".")
                mutationStr = splitMut[1]

                foundFiles += 1

                for line in log:
                    if not line.isspace():
                        meta = line.split(":")
                        if len(meta) > 1 and meta[0] == "ENERGY":
                            orig = meta[1].split()
                            data = []
                            data.append(pdbID)
                            data.append(mutationStr)
                            data.append(orig[0])
                            data.append(orig[12])
                            
                            writer.writerow(data)
                            numRows += 1
                log.close()

        out.close()
        print "Aggregated {0} logs into {1}, {2} rows".format(foundFiles, outFileName, numRows)