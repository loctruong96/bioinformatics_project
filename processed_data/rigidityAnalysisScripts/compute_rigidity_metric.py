import mpu.pd
import csv

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


baselines_files = ['base_2DM8_raw_ra_cluster_data.csv',
                   'base_3RL7_raw_ra_cluster_data.csv',
                   'base_MBDLG_raw_ra_cluster_data.csv']
mutation_files = ['2DM8_raw_ra_cluster_data.csv',
                  '3RL7_raw_ra_cluster_data.csv',
                  'MBDLG_raw_ra_cluster_data.csv']

baselines = process_csv(baselines_files)
mutation = process_csv(mutation_files)
results = {}
for PDB in baselines:
    WT = baselines[PDB]['G2837S']
    results[PDB] = {}
    for mut in mutation[PDB]:
        LRC = max(mutation[PDB][mut].keys())
        results[PDB][mut] = 0
        for i in range(1, LRC+1):
            WT_i = 0
            Mut_i = 0
            if i in mutation[PDB][mut]:
                Mut_i = mutation[PDB][mut][i]
            if i in WT:
                WT_i = WT[i]
            results[PDB][mut] += i * (WT_i - Mut_i)

csv_header = ['PDB', 'MUTATION', 'RIGIDITY_METRIC']
for PDB in results:
    csv_rows = []
    for mut in results[PDB]:
        row = [PDB, mut, results[PDB][mut]]
        csv_rows.append(row)
    with open(f'{PDB}_rigidity_metric.csv', 'w', newline='') as csv_file:
        wr = csv.writer(csv_file, quoting=csv.QUOTE_NONNUMERIC)
        wr.writerow(csv_header)
        for row in csv_rows:
            wr.writerow(row)
