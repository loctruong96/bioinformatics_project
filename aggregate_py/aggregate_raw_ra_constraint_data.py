import sys
import os
import csv


#Input configuration
dir = "./preliminary_exp_output/multiMutant_out/"
pdbID = sys.argv[1].upper()
tmp = os.listdir(dir)
mutName = ""
for folder in tmp:
    if folder[:len(pdbID)] == pdbID:
        mutName = folder
        break
dir+=mutName
#example dir: "/preliminary_exp_output/multiMutant_out/2DM8G2837:2843_out/"

#Output configuration
outFilename = pdbID+"_raw_ra_bond_data.csv"
columns = ["PDBID", "MUTATION", "CONSTRAINT", "COUNT"]

outDir = "processed_data/"+pdbID

if not os.path.exists("processed_data"):
    os.makedirs("processed_data")
if not os.path.exists(outDir):
    os.makedirs(outDir)


out = open(outDir+"/"+outFilename,"w+")
writer = csv.writer(out, delimiter=",")
writer.writerow(columns)

numRows = 0

mutations = os.listdir(dir)
for mutation in mutations:
    fileNameRoot = "{0}/{1}/{2}.all.ra.out_RA/{2}.all.processed.pdb/user/{2}.all.processed.pdb".format(dir,mutation, pdbID)
    if len(pdbID) > 4:
        fileNameRoot = "{0}/{1}/{3}.all.ra.out_RA/{3}.all.processed.pdb/user/{3}.all.processed.pdb".format(dir,mutation,pdbID,pdbID[:4])
    filename = fileNameRoot+"_MetricsPDB.txt"
    if os.path.exists(filename):
        log = open(filename, "r") 
        mutationStr = mutation.split(".")[1]
        reachedContstraitTypes = False
        for line in log:
            line = line[:-1].strip()
            if( not reachedContstraitTypes):
                if line == "=========================":
                    reachedContstraitTypes = True
            else:
                split = line.split(": ")
                data = [pdbID, mutationStr, split[0], split[1]]
                writer.writerow(data)
                numRows += 1
        log.close()

out.close()
print "Aggregated {0} logs into {1}, {2} rows".format(len(mutations), outFilename, numRows)