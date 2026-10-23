"""Read-only source analysis for the September 13 lab update. No model refit."""
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ECHO = ROOT / 'protocols/echo'
SNAP = ECHO / 'experiment_setup_output/20260909_paths/20260909T220918_350704Z'
sys.path.insert(0, str(ECHO))
from acquisition_round_20260909 import read_matrix, load_history


def icc(a):
    n, k = a.shape
    between = k * np.var(a.mean(axis=1), ddof=1)
    within = np.square(a-a.mean(axis=1)[:, None]).sum() / (n*(k-1))
    return float((between-within)/(between+(k-1)*within))


def describe(d, values):
    d = d.copy()
    d['raw'] = d['Destination Well'].map(values)
    reps = d.pivot(index='community_mask', columns='technical_replicate', values='raw')
    edge = [w for w in values.index if (w[0] in 'AP' or int(w[1:]) == 24) and int(w[1:]) != 1]
    pears, ranks = reps.corr(), reps.corr(method='spearman')
    meds = reps.median(axis=1)
    out = {
        'culture_median': float(d.raw.median()), 'edge_excluding_column1_median': float(values.loc[edge].median()),
        'plate_median': float(values.median()), 'icc': icc(reps.to_numpy()),
        'pearson_pairs': [float(pears.loc[a,b]) for a,b in [(1,2),(1,3),(2,3)]],
        'spearman_pairs': [float(ranks.loc[a,b]) for a,b in [(1,2),(1,3),(2,3)]],
        'median_replicate_range': float((reps.max(axis=1)-reps.min(axis=1)).median()),
        'largest_single_well': float(d.raw.max()), 'best_triplicate_median': float(meds.max()),
        'column1': [float(values[f'{r}1']) for r in 'ABCDEFGHIJKLMNOP'],
        'row_medians': [float(d.loc[d['Destination Well'].str.startswith(r), 'raw'].median()) for r in 'BCDEFGHIJKLMNO'],
        'replicates': reps.reset_index().rename(columns={1:'rep1',2:'rep2',3:'rep3'}).to_dict('records'),
        'community_medians': meds.tolist(),
    }
    return out, d


def main():
    design = pd.read_csv(SNAP/'12strains/design.csv')
    candidates = pd.read_csv(SNAP/'12strains/candidates.csv').set_index('community_mask')
    history, hist_wells, centers, files = load_history(ECHO, 'single')
    data = {'new': [], 'previous': [], 'history_counts': {'run_observations':len(history), 'unique':int(history.community_mask.nunique())}}
    files += [SNAP/'12strains/design.csv', SNAP/'12strains/candidates.csv', SNAP/'12strains/source_demand.csv', SNAP/'12strains/source_plan.csv']
    for path in sorted((ECHO/'readout/11-09').glob('*.csv')):
        values, meta = read_matrix(path)
        plate = int(re.search(r'ID1: plate_(\d+)', meta)[1])
        d = design[design['Destination Plate Name'].str.endswith(f'plate_{plate}')]
        assert len(d) == 306 and d.community_mask.nunique() == 102
        out, wells = describe(d, values)
        out.update({'plate':plate,'path':str(d['Source Strategy'].iloc[0]),'file':str(path.relative_to(ROOT))})
        med = wells.groupby('community_mask').raw.median()
        col = 'ei_posterior_mean' if plate == 1 else 'tm_posterior_mean'
        out['prediction_spearman'] = float(med.rank().corr(candidates.loc[med.index, col].rank()))
        top = wells.groupby('community_mask').agg(members=('members','first'), median=('raw','median'),
            minimum=('raw','min'),maximum=('raw','max'),wells=('Destination Well',lambda x:'|'.join(x)))
        out['top_descriptive_only'] = top.sort_values('median',ascending=False).head(5).reset_index().to_dict('records')
        data['new'].append(out)
        files.append(path)
    for plate, arm in enumerate(['bayesian_ei','max_diversity','bayesian_posterior_mean','random'],1):
        d = hist_wells[hist_wells.run_label.eq('20260902') & hist_wells.path.eq(arm)].copy()
        values, _ = read_matrix(ECHO/'readout/04-09'/d.reader_file.iloc[0])
        out, _ = describe(d, values)
        out.update({'plate':plate, 'path':arm})
        data['previous'].append(out)
    demand = pd.read_csv(SNAP/'12strains/source_demand.csv')
    source = pd.read_csv(SNAP/'12strains/source_plan.csv',dtype={'lab_marking':str})
    total = demand.groupby('Source Well').transfer_demand_uL.sum()
    loading = []
    for s in source.to_dict('records'):
        w = s['Source Well']
        loading.append({'species':s['lab_marking'],'well':w,'all_four_demand_ul':float(total[w]),
            'single_fill_theoretical_remaining_ul':50-float(total[w]),
            'worst_per_plate_remaining_ul':float(demand.loc[demand['Source Well'].eq(w),'expected_remaining_uL'].min())})
    data['loading'] = loading
    data['standards_concentrations'] = [12 / 2**r for r in range(16)]
    data['standard_mapping_source'] = 'User confirmed 12 mM at A1 and halving with each row, on all four plates.'
    data['actual_transfer_log_available'] = False
    data['refill_between_plates_confirmed'] = False
    assert sum(r['all_four_demand_ul'] > 50 for r in loading) == 7
    assert sum(r['all_four_demand_ul'] > 35 for r in loading) == 9
    screen = pd.read_csv(ROOT/'results/screening/genome_scores.tsv',sep='\t')
    ldh = pd.read_csv(ROOT/'results/screening/ldh_adjudication.tsv',sep='\t')
    data['screen'] = {'genomes':len(screen),'proteins_adjudicated':len(ldh),
        'genomes_by_call':{k:int(v) for k,v in ldh.groupby('call').genome.nunique().items()}}
    files += [ROOT/'results/screening/genome_scores.tsv', ROOT/'results/screening/ldh_adjudication.tsv']
    # Verify the same-plate comparator directly from the prior ordinary readout.
    prior = hist_wells[hist_wells.run_label.eq('20260902') & hist_wells.path.eq('max_diversity')]
    members = {'294','196','291','264'}
    community = prior[prior.members.str.split('|').map(lambda xs:set(xs)==members)]
    mono = prior[prior.members.isin(members)]
    data['candidate_comparison'] = {'members':sorted(members), 'community':community.raw_absorbance.tolist(),
        'monocultures':{s:mono.loc[mono.members.eq(s),'raw_absorbance'].tolist() for s in sorted(members)}}
    assert len(community)==3 and len(mono)==12
    data['source_hashes'] = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (HERE/'.build').mkdir(exist_ok=True)
    (HERE/'.build/data.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    lines = ['# Readout review: 11 September 2026','',
        'Analysis: 13 September. Original files and acquisition models remain unchanged.', '',
        '## Interpretation','',
        'Most culture-well readings are close to the edge reference. Replicate agreement is weak. '
        'Reported Echo assembly warnings mean intended composition is not verified. '
        'Do not treat a low reading as a confirmed low-producing community or train the next model on this round without resolving assembly status.', '',
        'The user confirmed standards in column 1: A1 = 12 mM, B1 = 6 mM, C1 = 3 mM, halving each row. '
        'All column-1 wells are excluded from the descriptive edge reference. '
        'The top standards respond on every plate, but lower levels overlap the background and spatial variation. '
        'No true zero-standard well was specified, and standard-versus-sample dilution equivalence is unconfirmed. '
        'No concentration conversion, arbitrary outlier deletion, or failed-well imputation is performed.', '',
        '| Plate | Culture median | Edge reference | ICC(1,1) | Median replicate range |',
        '|---|---:|---:|---:|---:|']
    for p in data['new']:
        lines.append(f"| {p['plate']} | {p['culture_median']:.4f} | {p['edge_excluding_column1_median']:.4f} | {p['icc']:.3f} | {p['median_replicate_range']:.4f} |")
    lines += ['', 'ICC compares between-community variation with within-triplicate variation. Negative values are retained.', '',
        '## Source loading audit','',
        'The historical statement that all wells retain at least 26 µL is correct **per destination plate with a new 50 µL fill**. '
        'It is not true for all four plates from one initial fill. Seven sources require more than 50 µL across all four plates. '
        'Nine exceed the 35 µL usable volume associated with the previous 15 µL reserve. '
        'The user confirmed there were no refills. The planning check therefore used the wrong loading assumption for the actual run. '
        'The cumulative demand explains why the run could not complete from one fill, even if all wells initially contained exactly 50 µL. '
        'No instrument report has been found in the scoped project files, so the exact failed transfers and affected compositions remain unresolved. '
        'This does not exclude additional assay or instrument effects.', '',
        '| Source | Species | Demand across four plates (µL) | Minimum remainder per plate after 50 µL fill |',
        '|---|---:|---:|---:|']
    for r in loading:
        lines.append(f"| {r['well']} | {r['species']} | {r['all_four_demand_ul']:.3f} | {r['worst_per_plate_remaining_ul']:.3f} |")
    lines += ['', '## Earlier evidence worth retaining','',
        f"The September 4 four-member candidate 294 + 196 + 291 + 264 has ordinary readings {data['candidate_comparison']['community']}. "
        f"The highest constituent-monoculture reading on the same plate is {mono.raw_absorbance.max():.3f}. "
        'This is candidate evidence, not proof of synergy: starting total biomass increases with richness and no independent validation is established.', '',
        '## Needed to resolve this round','',
        'Confirm blank identities, standard/sample assay dilutions and plate execution order. '
        'Match Echo survey/transfer exceptions to each requested transfer before deciding which compositions were assembled. '
        'Inspect source plate type and liquid-class settings against the actual labware. '
        'Do not infer successful assembly from the commanded picklist alone.', '', '## Sources','']
    lines += [f'- `{name}` (SHA-256 `{digest}`)' for name,digest in data['source_hashes'].items()]
    (HERE/'readout_review.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'new':[{k:p[k] for k in ['plate','culture_median','edge_excluding_column1_median','icc']} for p in data['new']],
        'candidate':data['candidate_comparison'],'screen':data['screen']},indent=2))


if __name__ == '__main__':
    main()
