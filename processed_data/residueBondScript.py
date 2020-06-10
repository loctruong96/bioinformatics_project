import os
import csv
import numpy as np

def getLigandAtoms(pdb, residueNumArray):
    # returns len(residueNumArray) size array of 1 or 0
    residues = {}
    for line in pdb:
        if(line[0:4] == "ATOM"):
            resNum = int(line[22:26])
            if(resNum >= residueNumArray[0] and resNum <= residueNumArray[1]):
                if resNum not in residues.keys():
                    residues[resNum] = []
                residues[resNum].append(int(line[6:11]))
    return residues  
                
def checkBonds(bonds, residueAtoms):
    isBonded = {}
    for residue in residueAtoms:
        isBonded[residue] = 0
    
    for line in bonds:
        bond = line.split()
        leftAtom = int(bond[1])
        rightAtom = int(bond[3])
        for residue in residueAtoms:
            if (leftAtom in residueAtoms[residue] and rightAtom not in residueAtoms[residue]) or (leftAtom not in residueAtoms[residue] and rightAtom in residueAtoms[residue]):
                isBonded[residue] = 1
                break
    return isBonded




groups = {
    "mbDLG_1": ["2DM8", "3RL7", "MBDLG"],
    "mbDLG_2": ["2BYG", "MBDLG2HB008", "MBDLG2HB010"],
    "mbDLG_3": ["2HE2", "MBDLG3"]
}

ligand = {
    "mbDLG_1": [2837, 2843],
    "mbDLG_2": [2838, 2843],
    "mbDLG_3": [513, 517]
}

ligandChain = {
    "mbDLG_1": "G",
    "mbDLG_2": "B",
    "mbDLG_3": "B"
}


out = open("ligandBonds.csv", "w+")
outwriter = csv.writer(out, delimiter=",")
outwriter.writerow(["PBDID", "MUTATION", "RESIDUE", "HBONDED", "HPHOBIC"])


for group in groups.keys():
    for pdbID in groups[group]:
        
        dir = "../raw_output/"+group+"_out/baselines/multiMutant_out/"

        mutName = ""
        tmp = os.listdir(dir)
        for folder in tmp:
            if folder[:len(pdbID)] == pdbID:
                mutName = folder
                break

        dir = dir+mutName+"/"
        tmp = os.listdir(dir)
        for folder in tmp:
            mutName = folder
            break

        dir = dir+mutName+"/"

        mutantDir = dir+pdbID[0:4]+".all.cur.out_C/"
        pdbFile = mutantDir+pdbID[0:4]+".all.processed.pdb.knr"
        #print(pdbFile)
        hbondFile = mutantDir+pdbID[0:4]+".all.HBonds.bnd.knr"
        hphobeFile = mutantDir+pdbID[0:4]+".all.HPhobes.bnd.knr"

        hbonds = open(hbondFile, "r")
        hphobes = open(hphobeFile, "r")
        pdb = open(pdbFile, "r")

        residueAtoms = getLigandAtoms(pdb, ligand[group])
        pdb.close()

        residueHbonded = checkBonds(hbonds, residueAtoms)
        hbonds.close()

        residueHPhobe = checkBonds(hphobes, residueAtoms)
        hphobes.close()


        for residue in residueAtoms:
            outwriter.writerow([pdbID, "WT", ligandChain[group]+str(residue), residueHbonded[residue],residueHPhobe[residue]])


        dir = "../raw_output/"+group+"_out/multiMutant_out/"
        
        tmp = os.listdir(dir)
        for folder in tmp:
            if folder[:len(pdbID)] == pdbID:
                mutName = folder
                break

        dir = dir+mutName+"/"



        mutants = os.listdir(dir)


        for mutant in mutants:
            if mutant != "sequentialPipelineInvocation.sh":
                #pdbFile = dir+mutant+"/"+mutant+"_em.pdb"
                #print(pdbFile)
                mutantDir = dir+mutant+"/"+pdbID[0:4]+".all.cur.out_C/"
                pdbFile = mutantDir+pdbID[0:4]+".all.processed.pdb.knr"
                hbondFile = mutantDir+pdbID[0:4]+".all.HBonds.bnd.knr"
                hphobeFile = mutantDir+pdbID[0:4]+".all.HPhobes.bnd.knr"

                hbonds = open(hbondFile, "r")
                hphobes = open(hphobeFile, "r")
                pdb = open(pdbFile, "r")

                mutation = mutant.split(".")[1]

                residueAtoms = getLigandAtoms(pdb, ligand[group])
                pdb.close()

                residueHbonded = checkBonds(hbonds, residueAtoms)
                hbonds.close()

                residueHPhobe = checkBonds(hphobes, residueAtoms)
                hphobes.close()


                for residue in residueAtoms:
                    outwriter.writerow([pdbID,mutation, ligandChain[group]+str(residue), residueHbonded[residue],residueHPhobe[residue]])


out.close()
