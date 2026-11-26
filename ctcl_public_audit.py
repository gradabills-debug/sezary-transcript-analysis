"""Focused audit of public processed CTCL data; not a reproduction of clustering."""
from pathlib import Path
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import io
import json
import hashlib
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'ctcl_public_data'
DATA.mkdir(exist_ok=True)
GEX = ['GSM8200498', 'GSM8200500', 'GSM8200502']
GENES = ['CD3D','CD4','CD8A','CD8B','PTGDR2','HPGDS','CRLF2','KLRB1']
VDJ = {'GSM8200498':'GSM8200499', 'GSM8200500':'GSM8200501', 'GSM8200502':'GSM8200503'}
HLA08_PAIR = (('TRA','CAVTGGGNKLTF'), ('TRB','CASSFWTGGNTEAFF'))

def fetch(url, path):
    if path.exists() and path.stat().st_size > 0:
        return path
    with urlopen(url, timeout=45) as response:
        tmp = path.with_suffix(path.suffix + '.part')
        with tmp.open('wb') as handle:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                handle.write(block)
        tmp.replace(path)
    return path

def get_sample(accession):
    meta = DATA / (accession + '.soft.txt')
    fetch('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + accession + '&targ=self&form=text&view=full', meta)
    lines = meta.read_text(encoding="utf-8").splitlines()
    title = next(x.split(' = ', 1)[1] for x in lines if x.startswith('!Sample_title'))
    urls = [x.split(' = ', 1)[1].replace('ftp://', 'https://') for x in lines if x.startswith('!Sample_supplementary_file')]
    wanted = [x for x in urls if any(x.endswith(s) for s in ('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz', 'filtered_contig_annotations.csv.gz'))]
    paths = {}
    for url in wanted:
        path = DATA / url.rsplit('/', 1)[1]
        fetch(url, path)
        key = 'contigs' if path.name.endswith('filtered_contig_annotations.csv.gz') else next(s for s in ('barcodes', 'features', 'matrix') if ('_' + s + '.') in path.name)
        paths[key] = str(path)
    result = {'accession': accession, 'title': title, 'files': paths}
    print(json.dumps({'downloaded': accession, 'title': title, 'bytes': sum(Path(p).stat().st_size for p in paths.values())}), flush=True)
    return result

def analyze(samples, vdj_samples):
    count_rows, expression_rows, clone_rows, cells_out = [], [], [], []
    for sample in samples:
        accession = sample['accession']
        features = pd.read_csv(sample['files']['features'], sep='\t', header=None)
        barcodes = pd.read_csv(sample['files']['barcodes'], sep='\t', header=None)[0].astype(str).to_numpy()
        with gzip.open(sample['files']['matrix'], 'rb') as handle, threadpool_limits(limits=2):
            matrix = mmread(handle).tocsc()
        matrix.eliminate_zeros()
        symbols = features[1].astype(str).to_numpy()
        total = np.asarray(matrix.sum(axis=0)).ravel()
        detected = np.diff(matrix.indptr)
        mt = np.asarray(matrix[np.char.startswith(symbols.astype(str), 'MT-'),:].sum(axis=0)).ravel()
        mito = np.divide(100 * mt, total, out=np.zeros_like(total, dtype=float), where=total>0)
        # Strict upper bound reproduces the submitted SS counts.
        # An inclusive 6000-gene cutoff retains one extra SS cell (HLA19).
        qc = (detected>=300) & (detected<6000) & (mito<=15)
        raw, expression = {}, {}
        for gene in GENES:
            row = np.flatnonzero(symbols==gene)
            if len(row) != 1:
                raise ValueError((accession,gene,len(row)))
            raw[gene] = np.asarray(matrix[row[0],:].toarray()).ravel()
            expression[gene] = np.log1p(np.divide(10000 * raw[gene], total, out=np.zeros_like(total,dtype=float), where=total>0))
        masks = {'qc':qc}
        for gene in ['CD3D']:
            masks[gene+'_raw_gt0'] = qc & (raw[gene]>0)
            masks[gene+'_raw_gt1'] = qc & (raw[gene]>1)
            masks[gene+'_lognorm_gt1'] = qc & (expression[gene]>1)
        selected = masks['CD3D_lognorm_gt1']
        count = {'gex_accession':accession, 'title':sample['title'], 'matrix_cells':len(barcodes)}
        count.update({k:int(v.sum()) for k,v in masks.items()})
        cell = pd.DataFrame({'gex_accession':accession,'barcode':barcodes,'qc_pass':qc,'selected_CD3D_lognorm_gt1':selected,'n_genes':detected,'n_counts':total,'pct_mito':mito})
        for gene in GENES:
            cell[gene+'_count'] = raw[gene]
            cell[gene+'_lognorm'] = expression[gene]
        paired = np.zeros(len(barcodes),dtype=bool)
        target = paired.copy()
        local_top = paired.copy()
        pairs, raw_ids = {}, {}
        if accession in VDJ:
            vdjpath = vdj_samples[VDJ[accession]]['files']['contigs']
            vdj = pd.read_csv(vdjpath)
            vdj = vdj[vdj['productive'] & vdj['high_confidence'] & vdj['is_cell'] & vdj['chain'].isin(['TRA','TRB'])]
            for barcode, group in vdj.groupby('barcode',sort=False):
                chains = tuple(sorted(set(zip(group['chain'],group['cdr3']))))
                if {'TRA','TRB'}.issubset({x[0] for x in chains}):
                    pairs[barcode] = chains
                    raw_ids[barcode] = ';'.join(sorted(set(group['raw_clonotype_id'].dropna())))
            selected_pairs = Counter(pairs[b] for b,s in zip(barcodes,selected) if s and b in pairs)
            top = selected_pairs.most_common(1)[0][0] if selected_pairs else None
            paired = np.array([b in pairs for b in barcodes])
            target = np.array([pairs.get(b)==HLA08_PAIR for b in barcodes])
            local_top = np.array([pairs.get(b)==top for b in barcodes])
            cell['paired_cdr3_aa'] = [json.dumps(pairs[b]) if b in pairs else '' for b in barcodes]
            cell['raw_clonotype_id_local'] = [raw_ids.get(b,'') for b in barcodes]
            cell['is_HLA08_target_pair'] = target
            cell['is_library_top_pair'] = local_top
            count['all_VDJ_paired_cells_before_GEX_QC'] = len(pairs)
            for maskname, mask in masks.items():
                count[maskname+'_paired'] = int((mask & paired).sum())
                count[maskname+'_HLA08_target'] = int((mask & target).sum())
            for pair, n in selected_pairs.most_common(10):
                clone_rows.append({'gex_accession':accession,'paired_cdr3_aa':json.dumps(pair),'selected_cells':n,'selected_T_cells':int(selected.sum()),'selected_paired_T_cells':int((selected&paired).sum()),'percent_all_selected_T':100*n/selected.sum(),'percent_paired_selected_T':100*n/(selected&paired).sum()})
        groups = {'all_selected':selected}
        if accession in VDJ:
            groups.update({'library_top_pair':selected&local_top,'other_T_cells':selected&~local_top,'other_paired_T_cells':selected&paired&~local_top,'other_CD4_detected_T_cells':selected&~local_top&(raw['CD4']>0),'HLA08_target_pair':selected&target})
        for groupname, mask in groups.items():
            n = int(mask.sum())
            if not n:
                continue
            for gene in GENES:
                x = expression[gene][mask]
                expression_rows.append({'gex_accession':accession,'group':groupname,'gene':gene,'n_cells':n,'n_detected':int((raw[gene][mask]>0).sum()),'pct_detected':100*float((raw[gene][mask]>0).mean()),'mean_lognorm':float(x.mean()),'median_lognorm':float(np.median(x))})
            for label, genes in [('PTGDR2_HPGDS_CRLF2',['PTGDR2','HPGDS','CRLF2']),('PTGDR2_HPGDS',['PTGDR2','HPGDS']),('HPGDS_CRLF2',['HPGDS','CRLF2'])]:
                coexp = np.logical_and.reduce([raw[g]>0 for g in genes])
                expression_rows.append({'gex_accession':accession,'group':groupname,'gene':'COEXPR_'+label,'n_cells':n,'n_detected':int((coexp&mask).sum()),'pct_detected':100*float(coexp[mask].mean()),'mean_lognorm':None,'median_lognorm':None})
        count_rows.append(count)
        cells_out.append(cell)
        print('ANALYZED '+json.dumps(count),flush=True)
        del matrix
    out = ROOT/'ctcl_audit_results'
    out.mkdir(exist_ok=True)
    pd.DataFrame(count_rows).to_csv(out/'sample_counts.csv',index=False)
    pd.DataFrame(expression_rows).to_csv(out/'marker_expression.csv',index=False)
    pd.DataFrame(clone_rows).to_csv(out/'top_clonotypes.csv',index=False)
    pd.concat(cells_out,ignore_index=True).to_csv(out/'cell_marker_audit.csv.gz',index=False,compression='gzip')
    print('ANALYSIS_COMPLETE',flush=True)

if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=4) as executor:
        samples = list(executor.map(get_sample, GEX))
    with ThreadPoolExecutor(max_workers=3) as executor:
        vdj_samples = {s['accession']:s for s in executor.map(get_sample, VDJ.values())}
    (DATA / 'samples.json').write_text(json.dumps(samples, indent=2))
    print('DOWNLOAD_COMPLETE', flush=True)
    analyze(samples, vdj_samples)
    manifest=[]
    for sample in samples + list(vdj_samples.values()):
        for key,path in sample['files'].items():
            p=Path(path)
            manifest.append({'accession':sample['accession'],'file_kind':key,'file_name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (ROOT/'ctcl_audit_results'/'source_manifest.json').write_text(json.dumps(manifest,indent=2))

    import sys
    versions = {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__,
                'scipy': __import__('scipy').__version__,
                'threadpoolctl': __import__('threadpoolctl').__version__}
    (ROOT/'ctcl_audit_results'/'runtime_versions.json').write_text(json.dumps(versions, indent=2), encoding='utf-8')
