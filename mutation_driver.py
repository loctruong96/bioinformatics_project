import os
from tqdm import tqdm
import argparse

"""
check the output folder of multiMutant for the following
* if the expected mutant PDB files exists
* if the expected mutant energy minimized PDB file exists
* if each mutant has the correct mutated residue
* if each mutant pdb and output fasta matched.
* if energy minimization is performed by comparing WT pdb and mutated_em pdb's non mutated atom
* if missing fasta are only the duplicated
"""


def get_args():
    parser = argparse.ArgumentParser(description='Basic tests for multiMutant output.')
    parser.add_argument('--output_dir', help='multiMutant output directory.', required=True)
    parser.add_argument('--protein', help='name of protein, e.g 1HHP.', required=True)
    return parser.parse_args()


# Returns the number of atoms in the non-mutated WT residues that have changed positions in the mutation
def em_test(path_to_wild_type_pdb, path_to_mutated_em_pdb, mutation_residue):
    wt = extract_atom_records(path_to_wild_type_pdb)
    mt = extract_atom_records(path_to_mutated_em_pdb)

    if wt[-1][1] != mt[-1][1]:
        print(f'ERROR, Wild Type and Mutant have different residue counts ( {len(wt)} vs {len(mt)} )')
        return -1

    # create a dictionary of atoms, key is (atom name, resNum), value is (x, y, z)
    # Note, the keys and values are strings, with potentially leading and trailing spaces according to PDB file format
    atom_dict = {}
    for atom in mt:
        if int(atom[1]) != mutation_residue:
            atom_dict[(atom[0], atom[1])] = (atom[2], atom[3], atom[4])
    count = 0
    for wt_atom in wt:
        if (wt_atom[0], wt_atom[1]) in atom_dict:
            mt_pos = atom_dict.get((wt_atom[0], wt_atom[1]))
            if wt_atom[2] != mt_pos[0] or wt_atom[3] != mt_pos[1] or wt_atom[4] != mt_pos[2]:
                count += 1
    return count


# Returns a list of atom records
# all atom records contain the following fields:
# Atom name, resSeq, x, y, z
# All fields are strings, with leading and trailing spaces according to the PDB file format
def extract_atom_records(path_to_pdb):
    f = open(path_to_pdb, "r")
    line = f.readline()
    records = []
    while line:
        line = f.readline()
        record_type = line[0:6]
        if record_type == "ATOM  ":
            new_record = [line[12:16], line[22:26], line[30:38], line[38:46], line[46:54]]
            records.append(new_record)
    f.close()
    return records


def read_proper_fasta(filepath):
    sequence = ''
    assert os.path.exists(filepath), f'{filepath} does not exists.'
    with open(filepath) as file_in:
        for i, line in enumerate(file_in):
            if i != 0:
                sequence += line.strip().upper()
    return sequence


def mutation_phase_1_test(output_dir, protein):
    assert os.path.exists(output_dir), f'{output_dir} directory does not exists.'
    if os.path.exists('pymol.log'):
        os.system('rm pymol.log')
    os.chdir(output_dir)
    result_dirs = next(os.walk('.'))[1]
    failed = []
    fasta_count = 0
    load_original_pdb = False
    for d in tqdm(result_dirs):
        fasta = f'{d}/{d}.fasta.txt'
        pdb = f'{d}/{d}.pdb'
        em_pdb = f'{d}/{d}_em.pdb'
        if not load_original_pdb:
            os.system(f'wget -q -O "{protein}.pdb" "https://files.rcsb.org/download/{protein}.pdb"')
            assert os.path.exists(f"{protein}.pdb")
            os.system(f'echo load {protein}.pdb > pymol_script.pml')
            os.system(f'echo save {protein}.fasta >> pymol_script.pml')
            os.system('pymol pymol_script.pml -qc >> pymol.log')
            assert os.path.exists(f"{protein}.fasta")
            load_original_pdb = True
        original_fasta = read_proper_fasta(f"{protein}.fasta")

        # count from 0
        residue = int(d[6:][:2]) - 1
        mutation = d[6:][2:]
        if os.path.exists(fasta):
            fasta_count += 1
            with open(fasta) as file_in:
                for line in file_in:
                    fasta_c = line.strip().upper()
            assert fasta_c[residue] == mutation, 'mismatched mutation'
            assert os.path.exists(pdb), 'mutated pdb does not exists.'
            os.system(f'echo load {pdb} > pymol_script.pml')
            os.system(f'echo save test.fasta >> pymol_script.pml')
            os.system('pymol pymol_script.pml -qc >> pymol.log')
            pymol_fasta = read_proper_fasta('test.fasta')
            assert len(pymol_fasta) == len(fasta_c), 'mismatch fasta size'
            for k, char in enumerate(pymol_fasta):
                assert char == fasta_c[k], 'mismatched pymol and pipline fasta'
            assert os.path.exists(em_pdb), 'mutated_em pdb does not exists.'
            assert em_test(pdb, em_pdb, residue) != 0, 'Energy minimization did not occur.'
        else:
            # check if the failed fasta is indeed a mutation from the same amino acid. e.g from A -> A again.
            # if not something went wrong
            if original_fasta[residue] != mutation:
                failed.append((original_fasta[residue], mutation, residue + 1))
    if os.path.exists('pymol_script.pml'):
        os.system('rm pymol_script.pml')
    if os.path.exists('test.fasta'):
        os.system('rm test.fasta')
    return failed


if __name__ == '__main__':
    args = get_args()
    failed_mutation = mutation_phase_1_test(args.output_dir, args.protein)
    assert len(failed_mutation) == 0, f'the {len(failed_mutation)} mutations failed {failed_mutation}'





