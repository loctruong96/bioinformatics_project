import mpu.pd
import csv
import os
import statistics
from scipy.special import softmax
import numpy as np

def process_csv(list_of_csv):
    results = {}
    for file in list_of_csv:
        df = mpu.pd.pd.read_csv(file, delimiter=',')
        dicts = df.to_dict('records')
        for record in dicts:
            if record['PDBID'] not in results:
                results[record['PDBID']] = {}
            if record['MUTATION'] not in results[record['PDBID']]:
                results[record['PDBID']][record['MUTATION']] = {record['CLUSTER_SIZE']: record['COUNT']}
            results[record['PDBID']][record['MUTATION']][record['CLUSTER_SIZE']] = record['COUNT']
    return results


def get_ra_csv(pdbs=None, path='../mbDLG_1/'):
    baselines_files = []
    mutation_files = []
    for pdb in pdbs:
        tmp_baseline_file = os.path.join(path, f'{pdb}/base_{pdb}_raw_ra_cluster_data.csv')
        assert os.path.exists(tmp_baseline_file), f'{tmp_baseline_file} does not exists'
        tmp_mutation_file = os.path.join(path, f'{pdb}/{pdb}_raw_ra_cluster_data.csv')
        assert os.path.exists(tmp_mutation_file), f'{tmp_mutation_file} does not exists'
        baselines_files.append(tmp_baseline_file)
        mutation_files.append(tmp_mutation_file)
    return baselines_files, mutation_files


def compute_ra_metric(baselines=None, mutation=None, baseline_id=None):
    final_results = {}
    abs_final_results = {}
    for PDB in baselines:
        WT = baselines[PDB][baseline_id]
        final_results[PDB] = {}
        abs_final_results[PDB] = {}
        for mut in mutation[PDB]:
            LRC = max(mutation[PDB][mut].keys())
            final_results[PDB][mut] = 0
            abs_final_results[PDB][mut] = 0
            for i in range(1, LRC + 1):
                WT_i = 0
                Mut_i = 0
                if i in mutation[PDB][mut]:
                    Mut_i = mutation[PDB][mut][i]
                if i in WT:
                    WT_i = WT[i]
                final_results[PDB][mut] += i * (WT_i - Mut_i)
                abs_final_results[PDB][mut] += i * abs(WT_i - Mut_i)
    return final_results, abs_final_results


def write_csv(header=None, data=None, csv_path=None):
    with open(csv_path, 'w', newline='') as csv_file:
        wr = csv.writer(csv_file, quoting=csv.QUOTE_NONNUMERIC)
        wr.writerow(header)
        for row in data:
            wr.writerow(row)


def generate_prediction(group_path=None, baseline_id=None, targets_pdb=None):
    print(f'Processing {group_path}')
    baselines_files, mutation_files = get_ra_csv(pdbs=targets_pdb, path=group_path)
    baselines = process_csv(baselines_files)
    mutation = process_csv(mutation_files)
    final_results, abs_final_results = compute_ra_metric(baselines=baselines, mutation=mutation,
                                                         baseline_id=baseline_id)

    csv_header = ('PDB', 'MUTATION', 'RIGIDITY_METRIC', 'NORMALIZED_RD_METRIC', 'Z_SCORE')
    print(f'Best prediction by norm_ra {csv_header}')
    for pdb in final_results:
        csv_rows = []
        all_abs_items = final_results[pdb].items()
        all_abs_values = [item[1] for item in all_abs_items]
        all_abs_keys = [item[0] for item in all_abs_items]
        abs_stdev = statistics.stdev(all_abs_values)
        abs_mean = statistics.mean(all_abs_values)
        abs_softmax = softmax(all_abs_values)
        best = ("", -1, -1, -1)
        for i, mut in enumerate(all_abs_keys):
            z_score = (final_results[pdb][mut] - abs_mean) / abs_stdev
            ra_norm = abs_softmax[i]
            row = (pdb, mut, final_results[pdb][mut], ra_norm, z_score)
            if ra_norm > best[3]:
                best = row
            csv_rows.append(row)
        print(f'Best prediction by norm_ra {best}')
        csv_path = os.path.join(group_path, f'{pdb}/{pdb}_rigidity_metric.csv')
        write_csv(header=csv_header, data=csv_rows, csv_path=csv_path)


# parameters
paths = ['../mbDLG_1/', '../mbDLG_2/', '../mbDLG_3/']
baseline_ids = ['G2837S', 'B2838Y', 'B513H']
targets_pdbs = [['2DM8', '3RL7', 'MBDLG'], ['2BYG', 'MBDLG2HB008', 'MBDLG2HB010'], ['2HE2', 'MBDLG3']]
for i, pt in enumerate(paths):
    group_path = pt
    baseline_id = baseline_ids[i]
    targets_pdb = targets_pdbs[i]
    generate_prediction(group_path=group_path, baseline_id=baseline_id, targets_pdb=targets_pdb)
