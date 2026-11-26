"""Rebuild the SS submission summaries from the supplied cell-level export.

Python 3.10+ standard library only. No hypothesis testing or malignant-cell
classification is performed. UMI rarefaction is analytical, not a simulation.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import platform
import statistics
from collections import Counter

ROOT = Path(__file__).resolve().parent
SAMPLES = [('GSM8200500', 'HLA08 reaction 1'),
           ('GSM8200502', 'HLA08 reaction 2'), ('GSM8200498', 'HLA19')]
GENES = ['PTGDR2', 'HPGDS', 'CRLF2', 'KLRB1']
JOINT = ['PTGDR2', 'HPGDS', 'CRLF2']
TARGET = (('TRA', 'CAVTGGGNKLTF'), ('TRB', 'CASSFWTGGNTEAFF'))
DEPTH = 1000


def truth(value):
    return str(value).strip().lower() == 'true'


def sequence_key(value):
    return tuple(sorted(tuple(x) for x in json.loads(value))) if value else ()


def read_cells(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    seen = set()
    for row in rows:
        identity = row['gex_accession'], row['barcode']
        if identity in seen:
            raise ValueError(f'Duplicate accession and barcode: {identity}')
        seen.add(identity)
        row['selected'] = truth(row['selected_CD3D_lognorm_gt1'])
        row['key'] = sequence_key(row.get('paired_cdr3_aa', ''))
        row['target'] = row['key'] == TARGET
        row['N'] = int(float(row['n_counts']))
        for gene in GENES + ['CD4', 'CD8A', 'CD8B']:
            row[gene] = int(float(row[gene + '_count']))
    return rows


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def absence_probability(total, marker, depth):
    """C(total-marker, depth) / C(total, depth), evaluated stably."""
    if not (0 <= marker <= total and 0 <= depth <= total):
        raise ValueError((total, marker, depth))
    if marker == 0:
        return 1.0
    if total - marker < depth:
        return 0.0
    return math.exp(math.fsum(math.log1p(-depth / (total - j))
                              for j in range(marker)))


def joint_probability(total, counts, depth):
    """Probability of detecting all distinct genes after UMI subsampling."""
    if sum(counts) > total:
        raise ValueError('Marker UMIs exceed total UMIs')
    terms = [1.0]
    for size in range(1, len(counts) + 1):
        terms.extend((-1) ** size * absence_probability(total, sum(x), depth)
                     for x in itertools.combinations(counts, size))
    return max(0.0, min(1.0, math.fsum(terms)))


def validate_probability_functions():
    # Small exact combinatorial examples, independent of the numerical formula.
    for total in (6, 10, 20):
        for depth in range(1, total + 1):
            for marker in range(total + 1):
                exact = (math.comb(total - marker, depth) / math.comb(total, depth)
                         if total - marker >= depth else 0.0)
                assert abs(absence_probability(total, marker, depth) - exact) < 1e-12
    labels = [0, 0, 1, 2, 2, -1, -1, -1, -1]
    successes = sum({labels[i] for i in draw} >= {0, 1, 2}
                    for draw in itertools.combinations(range(9), 4))
    assert abs(joint_probability(9, [2, 1, 2], 4)
               - successes / math.comb(9, 4)) < 1e-12


def marker_summary(label, rows):
    n = len(rows)
    result = {'population': label, 'n_cells': n}
    for gene in GENES:
        positive = sum(row[gene] > 0 for row in rows)
        result[gene + '_positive'] = positive
        result[gene + '_pct'] = 100 * positive / n
    joint = sum(all(row[gene] > 0 for gene in JOINT) for row in rows)
    result.update(joint_positive=joint, joint_pct=100 * joint / n)
    return result


def sensitivity_summary(label, rows, subset, depth):
    if not rows:
        raise ValueError(f'Empty sensitivity population: {label}')
    if min(row['N'] for row in rows) < depth:
        raise ValueError('Rarefaction depth exceeds at least one cell total')
    result = {'subset': subset, 'population': label, 'n_cells': len(rows),
              'median_total_UMIs': statistics.median(row['N'] for row in rows),
              'rarefaction_UMIs_per_cell': depth}
    for gene in GENES:
        result[gene + '_observed_positive'] = sum(row[gene] > 0 for row in rows)
        result[gene + '_observed_pct'] = 100 * result[gene + '_observed_positive'] / len(rows)
        result[gene + '_expected_pct'] = 100 * math.fsum(
            1 - absence_probability(row['N'], row[gene], depth) for row in rows) / len(rows)
    result['joint_observed_positive'] = sum(all(row[g] > 0 for g in JOINT) for row in rows)
    result['joint_observed_pct'] = 100 * result['joint_observed_positive'] / len(rows)
    result['joint_expected_pct'] = 100 * math.fsum(
        joint_probability(row['N'], [row[g] for g in JOINT], depth)
        for row in rows) / len(rows)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells', type=Path, default=ROOT / 'data/cell_marker_audit.csv.gz')
    parser.add_argument('--out', type=Path, default=ROOT / 'results')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    validate_probability_functions()
    cells = read_cells(args.cells)
    selected = [r for r in cells if r['selected'] and r['gex_accession']
                in {x[0] for x in SAMPLES}]
    samples, figure1, clones = [], [], []
    by_sample = {}
    for accession, label in SAMPLES:
        rows = [r for r in selected if r['gex_accession'] == accession]
        paired = [r for r in rows if r['key']]
        counts = Counter(r['key'] for r in paired)
        leading = counts.most_common(2 if accession == 'GSM8200498' else 1)
        by_sample[accession] = rows
        samples.append({'gex_accession': accession, 'population': label,
                        'selected_cells': len(rows), 'paired_cells': len(paired),
                        'paired_pct': 100 * len(paired) / len(rows),
                        'leading_clonotype_cells': '; '.join(str(n) for _, n in leading)})
        common = dict(gex_accession=accession, selected_cells=len(rows),
                      paired_cells=len(paired), unpaired_cells=len(rows)-len(paired))
        figure1.append(dict(panel='a', label=label, **common, clonotype_cells='',
                            percentage=100 * len(paired) / len(rows), numerator=len(paired),
                            denominator=len(rows), paired_cdr3_aa=''))
        for rank, (key, n) in enumerate(leading, 1):
            clone_label = f'HLA19 pair {rank}' if accession == 'GSM8200498' else label
            figure1.append(dict(panel='b', label=clone_label, **common, clonotype_cells=n,
                                percentage=100*n/len(paired), numerator=n, denominator=len(paired),
                                paired_cdr3_aa=json.dumps(key)))
            clones.append({'gex_accession': accession, 'label': clone_label,
                           'paired_cdr3_aa': json.dumps(key), 'n_cells': n})

    h08 = [r for r in selected if r['gex_accession'] in {'GSM8200500', 'GSM8200502'}]
    h19 = by_sample['GSM8200498']
    h19_keys = Counter(r['key'] for r in h19 if r['key']).most_common(2)
    groups = [
        ('HLA08 clonotype reaction 1', 'HLA08 clonotype\nreaction 1',
         [r for r in h08 if r['gex_accession']=='GSM8200500' and r['target']]),
        ('HLA08 clonotype reaction 2', 'HLA08 clonotype\nreaction 2',
         [r for r in h08 if r['gex_accession']=='GSM8200502' and r['target']]),
        ('Other HLA08 paired-TCR cells', 'Other HLA08 cells\nwith paired TCRs',
         [r for r in h08 if r['key'] and not r['target']]),
        ('HLA19 pair 1', 'HLA19 pair 1', [r for r in h19 if r['key']==h19_keys[0][0]]),
        ('HLA19 pair 2', 'HLA19 pair 2', [r for r in h19 if r['key']==h19_keys[1][0]])]
    markers, figure2 = [], []
    for label, plot_label, rows in groups:
        result = marker_summary(label, rows)
        markers.append(result)
        key = json.dumps(rows[0]['key']) if label != 'Other HLA08 paired-TCR cells' else 'Multiple paired sequence keys'
        figure2.append(dict(population=label, plot_label=plot_label, n_cells=len(rows),
                            paired_cdr3_aa=key, **{k:v for k,v in result.items()
                                                if k not in {'population','n_cells'}}))
    assert [m['n_cells'] for m in markers] == [773,1311,595,238,236]
    assert [m['joint_positive'] for m in markers] == [118,179,0,1,1]
    assert len(selected)==7929 and len(h08)==4903
    paired_total = sum(bool(r['key']) for r in selected)
    h08_paired = sum(bool(r['key']) for r in h08)
    target_n = sum(r['target'] for r in h08)
    assert (paired_total,h08_paired,target_n)==(4657,2679,2084)

    sensitivity = []
    for subset in ['all_paired', 'CD4_detected_CD8A_CD8B_not_detected']:
        for accession, reaction in SAMPLES[:2]:
            rows = [r for r in h08 if r['gex_accession']==accession and r['key']]
            if subset != 'all_paired':
                rows = [r for r in rows if r['CD4']>0 and r['CD8A']==0 and r['CD8B']==0]
            for target, label in [(True, f'{reaction} target clonotype'),
                                  (False, f'{reaction} other paired-TCR cells')]:
                group = [r for r in rows if r['target']==target]
                sensitivity.append(sensitivity_summary(label, group, subset, DEPTH))
    summary = {
        'SS_selected_cells':len(selected), 'SS_paired_cells':paired_total,
        'HLA08_selected_cells':len(h08), 'HLA08_paired_cells':h08_paired,
        'HLA08_target_cells':target_n, 'target_pct_all_HLA08':100*target_n/len(h08),
        'target_pct_paired_HLA08':100*target_n/h08_paired,
        'rarefaction_UMIs_per_cell':DEPTH,
        'minimum_HLA08_paired_cell_UMIs':min(r['N'] for r in h08 if r['key']),
        'interpretation':'Descriptive; no independent malignant-cell verification or patient-level replication.',
        'input_file_sha256':hashlib.sha256(args.cells.read_bytes()).hexdigest()}
    write_csv(args.out/'Supplementary_Table_S1.csv', samples)
    write_csv(args.out/'Supplementary_Table_S2.csv', markers)
    write_csv(args.out/'Figure_1_source_data.csv', figure1)
    write_csv(args.out/'Figure_2_source_data.csv', figure2)
    write_csv(args.out/'clonotype_sequences.csv', clones)
    write_csv(args.out/'Supplementary_Table_S3.csv', [r for r in sensitivity if r['subset']=='all_paired'])
    write_csv(args.out/'Supplementary_Table_S4.csv', [r for r in sensitivity if r['subset']!='all_paired'])
    write_csv(args.out/'sensitivity_results.csv', sensitivity)
    (args.out/'verified_totals.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (args.out/'summary_runtime.json').write_text(json.dumps(
        {'python':platform.python_version(),'implementation':platform.python_implementation(),
         'summary_dependencies':'Python standard library only'},indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'summary':summary,'sensitivity':sensitivity},indent=2))


if __name__ == '__main__':
    main()
