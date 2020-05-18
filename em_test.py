
import os

# Returns the number of atoms in the non-mutated WT residues that have changed positions in the mutation
def em_test(path_to_wild_type_pdb, path_to_mutated_em_pdb, mutation_residue):
    wt = extractAtomRecords(path_to_wild_type_pdb)
    mt = extractAtomRecords(path_to_mutated_em_pdb)
    
    if(wt[-1][1] != mt[-1][1]):
        print(f'ERROR, Wild Type and Mutant have different residue counts ( {len(wt)} vs {len(mt)} )')
        return -1
    
    # create a dictionary of atoms, key is (atom name, resNum), value is (x, y, z)
    # Note, the keys and values are strings, with potentially leading and trailing spaces according to PDB file format
    atom_dict = {}
    for atom in mt:
        if(int(atom[1]) != mutation_residue):
            atom_dict[(atom[0],atom[1])] = (atom[2], atom[3], atom[4])
    
    count = 0
    #matches = 0
    
    for wt_atom in wt:
        if( (wt_atom[0],wt_atom[1]) in atom_dict):
            mt_pos = atom_dict.get((wt_atom[0],wt_atom[1]))
            #matches += 1
            if( wt_atom[2] != mt_pos[0] or\
                wt_atom[3] != mt_pos[1] or\
                wt_atom[4] != mt_pos[2]):
                count += 1
                
    #print(f'Atoms in WT: {len(wt)} \nAtoms in mutant: {len(mt)} \nAtom comparisons: {matches} \nAtoms moved: {count}')
    return count


# Returns a list of atom records
# all atom records contain the following fields:
# Atom name, resSeq, x, y, z
# All fields are strings, with leading and trailing spaces according to the PDB file format
def extractAtomRecords(path_to_pdb):
    f = open(path_to_pdb, "r")
    line = f.readline()
    records = []
    
    while(line):
        line = f.readline()
        recordType = line[0:6]
        if(recordType == "ATOM  "):
            newRecord = []
            newRecord.append(line[12:16])   # Atom name.
            newRecord.append(line[22:26])   # Residue sequence number.
            newRecord.append(line[30:38])   # Orthogonal coordinates for X in Angstroms.
            newRecord.append(line[38:46])   # Orthogonal coordinates for Y in Angstroms.
            newRecord.append(line[46:54])   # Orthogonal coordinates for Z in Angstroms.
            
            records.append(newRecord)
            
    f.close()
    return records   

print(em_test("./1HHP_A34F/1HHP.pdb", "./1HHP_A34F/1HHP.A34F.pdb", 34))