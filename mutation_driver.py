import os
from tqdm import tqdm
import argparse

"""
check the output folder of multiMutant for the following
* if the expected mutant PDB files exists
* if the expected mutant energy minimized PDB file exists
* if each mutant has the correct mutated residue
* if there exists only one mutation per fasta
* if each mutant pdb and output fasta matched
* if energy minimization is performed by comparing WT pdb and mutated_em pdb's non mutated atom
* if missing fasta are only the duplicated
"""
ALL_ACIDS = ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")


def get_args():
    parser = argparse.ArgumentParser(description='Basic tests for multiMutant output.')
    parser.add_argument('--output_dir', help='multiMutant output directory.', required=True)
    parser.add_argument('--protein', help='name of protein, e.g 1HHP.', required=True)
    parser.add_argument('--chain', help='chain of protein, e.g A.', required=True)
    parser.add_argument('--start_residue', type=int, help='start residue', required=True)
    parser.add_argument('--last_residue', type=int, help='last residue', required=True)
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
    termination = False
    while line and not termination:
        line = f.readline()
        record_type = line[0:6]
        if record_type == "ATOM  ":
            new_record = [line[12:16], line[22:26], line[30:38], line[38:46], line[46:54]]
            records.append(new_record)
        elif record_type == "TER   ":
            termination = True
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


def mutation_and_rigidity_test(output_dir=None, protein=None, chain=None, start=1, end=0, test_em=True,
                               test_rigidity=False):
    assert start >= 1, 'count from 0, start must >= 1'
    assert os.path.exists(output_dir), f'{output_dir} directory does not exists.'
    if os.path.exists('pymol.log'):
        os.system('rm pymol.log')
    os.chdir(output_dir)
    result_dirs = []
    mutations = []
    for i in range(start, end+1):
        for acid in ALL_ACIDS:
            output_path = f'{protein}.{chain}{i}{acid}'
            assert os.path.exists(output_path), f'{output_path} directory does not exists.'
            result_dirs.append(output_path)
            # (residue location, mutation)
            mutations.append((i, acid))
    failed = []
    fasta_count = 0
    load_original_pdb = False
    print(f'Begin multiMutant output test for {output_dir} with energy minimization: {test_em} '
          f'and rigidity analysis: {test_rigidity}')
    for z, d in enumerate(tqdm(result_dirs)):
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
        residue = mutations[z][0] - 1
        mutation = mutations[z][1]
        if os.path.exists(fasta):
            fasta_count += 1
            with open(fasta) as file_in:
                for line in file_in:
                    fasta_c = line.strip().upper()
            assert fasta_c[residue] == mutation, 'mismatched mutation'
            assert len(original_fasta) - len([item for i, item in enumerate(fasta_c) if i != residue
                                              and item == original_fasta[i]]) == 1, f'detected more than one mutation.'
            assert os.path.exists(pdb), 'mutated pdb does not exists.'
            os.system(f'echo load {pdb} > pymol_script.pml')
            os.system(f'echo save test.fasta >> pymol_script.pml')
            os.system('pymol pymol_script.pml -qc >> pymol.log')
            pymol_fasta = read_proper_fasta('test.fasta')
            assert len(pymol_fasta) == len(fasta_c), 'mismatch fasta size'
            for k, char in enumerate(pymol_fasta):
                assert char == fasta_c[k], 'mismatched pymol and pipline fasta'
            if test_em:
                assert os.path.exists(em_pdb), 'mutated_em pdb does not exists.'
                assert em_test(pdb, em_pdb, residue) != 0, 'Energy minimization did not occur.'
            if test_rigidity:
                metrics_bbh = f'{d}/{protein}.{chain}.ra.out_RA/{protein}.{chain}.processed.pdb/user/' \
                              f'{protein}.{chain}.processed.pdb_postPG_MetricsBBH.txt'
                metrics_pdb = f'{d}/{protein}.{chain}.ra.out_RA/{protein}.{chain}.processed.pdb/user/' \
                              f'{protein}.{chain}.processed.pdb_MetricsPDB.txt'
                assert os.path.exists(metrics_bbh), f'{metrics_bbh} does not exists.'
                assert os.path.exists(metrics_pdb), f'{metrics_pdb} metrics_pdb does not exists.'

        else:
            # check if the failed fasta is indeed a mutation from the same amino acid. e.g from A -> A again.
            # if not something went wrong
            if original_fasta[residue] != mutation:
                failed.append((original_fasta[residue], mutation, residue + 1))
    if os.path.exists('pymol_script.pml'):
        os.system('rm pymol_script.pml')
    if os.path.exists('test.fasta'):
        os.system('rm test.fasta')
    return failed, result_dirs


if __name__ == '__main__':
    args = get_args()
    failed_mutations, successful_mutation_dirs = mutation_and_rigidity_test(output_dir=args.output_dir,
                                                                            protein=args.protein.upper(),
                                                                            chain=args.chain.upper(),
                                                                            start=args.start_residue,
                                                                            end=args.last_residue,
                                                                            test_em=True,
                                                                            test_rigidity=True)
    assert len(failed_mutations) == 0, f'the {len(failed_mutations)} mutations failed {failed_mutations}'
    print(f'Found {len(successful_mutation_dirs)} successful mutations and {len(failed_mutations)} failed mutations.')





