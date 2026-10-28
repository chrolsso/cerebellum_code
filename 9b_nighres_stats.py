import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

#input files and settings
long_file = 'output/nighres_long.csv'      # from 9a
qc_file = 'output/nighres_qc.csv'          # from 9a
merged = 'output/clinical_merged.csv'      # from 8a
outdir = 'output/stats/'
exclude = []                               # ids to leave out
formula = 'y ~ age + patient'              # add + C(sex) if wanted
min_n = 8                                  # regions with fewer subjects are skipped
qc_mask_range = (100, 180)                 # cm3, only used to print a list of subjects to check
#%%

os.makedirs(outdir, exist_ok=True)
df = pd.read_csv(long_file)
df = df[~df['id'].isin(exclude)].copy()
df['patient'] = (df['status'] == 'Patient').astype(int)

# volumes as % of TIV, like the FastSurfer analysis
df['y'] = df['value']
vol = df['measure'].str.endswith('_volume_cm3')
df.loc[vol, 'y'] = 100 * df.loc[vol, 'value'] / df.loc[vol, 'tiv']
df.loc[vol, 'measure'] = df.loc[vol, 'measure'].str.replace('_cm3', '_tiv_pct')

#%% patient versus control per region and measure, adjusted for age
out = []
for (measure, region), d in df.groupby(['measure', 'region']):
    d = d[['y', 'age', 'sex', 'patient']].dropna(subset=['y', 'age'])
    if len(d) < min_n or d['patient'].nunique() < 2: continue
    fit = smf.ols(formula, data=d).fit()
    ci = fit.conf_int().loc['patient']
    control_mean = d.loc[d['patient'] == 0, 'y'].mean()
    out.append({'measure': measure, 'region': region, 'n': int(fit.nobs), 'control_mean': control_mean,
                'diff': fit.params['patient'], 'ci_low': ci[0], 'ci_high': ci[1], 'rel_diff_pct': 100 * fit.params['patient'] / control_mean,
                'p': fit.pvalues['patient'], 'p_age': fit.pvalues['age']})
res = pd.DataFrame(out)
res['p_fdr'] = res.groupby('measure')['p'].transform(lambda p: multipletests(p, method='fdr_bh')[1])   # over the regions of each measure
res.to_csv(outdir + 'nighres_patient_vs_control.csv', index=False)
print(res[res['region'].isin(['all', 'left', 'right'])].round(4).to_string(index=False))

lob = res[res['region'].str.startswith('lobule_')]
for measure in ['cortex_volume_tiv_pct', 'thickness_median_mm']:
    print('\n', measure, 'per lobule, sorted by p')
    print(lob[lob['measure'] == measure].sort_values('p')[['region', 'n', 'control_mean', 'diff', 'rel_diff_pct', 'p', 'p_fdr']].round(4).to_string(index=False))

#%% nighres against FastSurfer, same subjects
cb = df[(df['region'] == 'all') & (df['measure'] == 'cortex_volume_cm3')][['id', 'status', 'value']].rename(columns={'value': 'nighres_cortex_cm3'})
fs = pd.read_csv(merged)[['id', 'cerebellum_cortex']].rename(columns={'cerebellum_cortex': 'fastsurfer_cortex_cm3'})
both = cb.merge(fs, on='id').dropna()
both.to_csv(outdir + 'nighres_vs_fastsurfer.csv', index=False)
if len(both) > 2:
    r, p = stats.pearsonr(both['nighres_cortex_cm3'], both['fastsurfer_cortex_cm3'])
    print(f'\nnighres vs FastSurfer cortex volume: n={len(both)}, r={r:.2f}, p={p:.3g}')

#%% subjects to look at, not removed
qc = pd.read_csv(qc_file)
bad = qc[(qc['cb_mask_cm3'] < qc_mask_range[0]) | (qc['cb_mask_cm3'] > qc_mask_range[1]) | (qc['wm_object_components'] != 1)]
print('\nQC to check:', bad[['id', 'cb_mask_cm3', 'wm_object_components']].to_string(index=False) if len(bad) else 'none')
print('written to', outdir)
