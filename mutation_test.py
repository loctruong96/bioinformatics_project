import os
from tqdm import tqdm
# output_dir = "1HHPA34:45_out"
# if os.path.exists(output_dir):
#     os.system(f'rm -r {output_dir}')
# os.system('./multiMutant.sh 1HHP A 34:45 em')
# assert os.path.exists(output_dir), 'missing output folder'
result_dirs = next(os.walk('.'))[1]
failed = []
print(result_dirs)
count = 0
for d in tqdm(result_dirs):
    fasta = f'{d}/{d}.fasta.txt'
    pdb = f'{d}/{d}.pdb'
    em_pdb = f'{d}/{d}_em.pdb'
    # count from 0
    residue = int(d[6:][:2]) - 1
    mutation = d[6:][2:]

    if os.path.exists(fasta):
        count += 1
        with open(fasta) as file_in:
            for line in file_in:
                fasta_c = line.strip().upper()
        assert fasta_c[residue] == mutation, 'mismatched mutation'
        if os.path.exists('pymol_script.pml'):
            os.system('rm pymol_script.pml')
        assert os.path.exists(pdb), 'mutated pdb does not exists.'
        assert os.path.exists(em_pdb), 'mutated_em pdb does not exists.'
        os.system(f'echo load {pdb} >> pymol_script.pml')
        os.system(f'echo save test.fasta >> pymol_script.pml')
        os.system('pymol pymol_script.pml -qc')
        pymol_fasta = ''
        with open('test.fasta') as file_in:
            for i, line in enumerate(file_in):
                if i != 0:
                    pymol_fasta += line.strip().upper()
        assert len(pymol_fasta) == len(fasta_c), 'mismatch fasta size'
        for k, char in enumerate(pymol_fasta):
            assert char == fasta_c[k], 'mismatched pymol and pipline fasta'

    else:
        failed.append(fasta)

print(f'fasta count {count}, total {len(result_dirs)}.')
print(failed)