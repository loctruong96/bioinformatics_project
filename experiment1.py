import os
import time

start_index = 2837
end_index = 2843
chain = 'G'
sequences = ['2DM8', 'MBDLG', '3RL7']
start = time.time()
total = start
for sequence in sequences:
    print(f'processing {sequence}')
    output_folder_name = f'{sequence}{chain}{start_index}:{end_index}_out/'
    assert os.path.exists('multiMutant.sh'), 'multiMutant script is missing.'
    if os.path.exists(output_folder_name):
        print(f'removing previous results folder {output_folder_name}')
        os.system(f'rm -r {output_folder_name}')
    print(f'execute: ./multiMutant.sh {sequence} {chain} {start_index}:{end_index} -em')
    os.system(f'./multiMutant.sh {sequence} {chain} {start_index}:{end_index} -em')
    print(f'command took {round(time.time()-start)} seconds.')
    start = time.time()
    assert os.path.exists(output_folder_name), f'output folder for {sequence} is missing.'
    os.chdir(output_folder_name)
    assert os.path.exists('sequentialPipelineInvocation.sh'), 'sequentialPipelineInvocation script is missing.'
    print('execute: ./sequentialPipelineInvocation.sh')
    os.system(f'./sequentialPipelineInvocation.sh')
    print(f'command took {round(time.time() - start)} seconds.')
    start = time.time()
    os.chdir('..')

print(f'total time taken {round(time.time()-total)} seconds.')
