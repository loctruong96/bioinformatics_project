import os
import time

from output_test import read_proper_fasta, mutation_and_rigidity_test
mutation_residues = 'SYLVTSV'
chain = 'A'
sequences = ['2DM8', 'MBDLG', '3RL7']
start = time.time()
for sequence in sequences:
    print(f'processing {sequence}')
    fasta_file = f'promute/{sequence}.fasta'
    assert os.path.exists(fasta_file), f'{fasta_file} does not exists.'
    fasta_file_content = read_proper_fasta(fasta_file)
    assert mutation_residues in fasta_file_content, f'{mutation_residues} does not exists in {sequence}.'
    start_index = fasta_file_content.index(mutation_residues)
    end_index = fasta_file_content.index(mutation_residues)+(len(mutation_residues)-1)
    assert start_index <= len(fasta_file_content)-1
    assert end_index <= len(fasta_file_content) - 1
    test = ""
    for i in range(start_index, end_index+1):
        test += fasta_file_content[i]
    assert test == mutation_residues, 'start index and end index mismatched.'
    output_folder_name = f'{sequence}{chain}{start_index}:{end_index}_out/'
    assert os.path.exists('multiMutant.sh'), 'multiMutant script is missing.'
    if os.path.exists(output_folder_name):
        print(f'removing previous results folder {output_folder_name}')
        os.system(f'rm -r {output_folder_name}')
    print(f'execute: ./multiMutant.sh {sequence} {chain} {start_index}:{end_index} -em')
    os.system(f'./multiMutant.sh {sequence} {chain} {start_index}:{end_index} -em')
    print(f'command took {round(time.time()-start)} seconds.')
    assert os.path.exists(output_folder_name), f'output folder for {sequence} is missing.'
    os.chdir(output_folder_name)
    assert os.path.exists('sequentialPipelineInvocation.sh'), 'sequentialPipelineInvocation script is missing.'
    print('execute: ./sequentialPipelineInvocation.sh')
    os.system(f'./sequentialPipelineInvocation.sh')
    print(f'command took {round(time.time() - start)} seconds.')
    os.chdir('..')

print(f'total time taken {round(time.time()-start)} seconds.')
