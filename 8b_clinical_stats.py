import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

#input files and settings
merged = 'output/clinical_merged.csv' # from 8a
outdir = 'output/stats/'
control = 'Control'
min_n = 3 # groups with fewer subjects are left out of the models

# One formula per model. Column names are the lower-case versions of your Excel columns. Remove or add depending on analysis
models = {
    'cerebellum_by_group': 'cerebellum_total ~ C(group, Treatment("Control"))',
    'cerebellum_by_group_brain': 'cerebellum_total ~ C(group, Treatment("Control")) + brain_total',
    'brain_by_group':    'brain_total ~ C(group, Treatment("Control"))',
    # 'cerebellum_full': 'cerebellum_total ~ C(group, Treatment("Control")) + brain_total + age + sex',
}
correlate = [('cerebellum_total', 'brain_total')]
#%%

os.makedirs(outdir, exist_ok=True)
df = pd.read_csv(merged)
n = df.groupby('group')['id'].transform('size')
print(df.groupby('group').size().to_string())
small = sorted(df.loc[n < min_n, 'group'].unique())
if small: print('left out (fewer than', min_n, 'subjects):', small)
df = df[n >= min_n]

#%% group means
cols = sorted({c for f in models.values() for c in [f.split('~')[0].strip()]})
df.groupby('group')[cols].agg(['count', 'mean', 'std']).to_csv(outdir + 'group_means.csv')

#%% correlations, all subjects and per group
corr = []
for a, b in correlate:
    for g, d in [('all', df)] + list(df.groupby('group')):
        d = d[[a, b]].dropna()
        if len(d) < 3: continue
        r, p = stats.pearsonr(d[a], d[b]); rho, ps = stats.spearmanr(d[a], d[b])
        corr.append({'x': a, 'y': b, 'group': g, 'n': len(d), 'pearson_r': r, 'p': p, 'spearman_rho': rho, 'p_spearman': ps})
pd.DataFrame(corr).to_csv(outdir + 'correlations.csv', index=False)

# regression + ANOVA on the same fitted model
coef_all = []
for name, formula in models.items():
    fit = smf.ols(formula, data=df).fit()
    anova = sm.stats.anova_lm(fit, typ=2)
    anova.to_csv(outdir + f'anova_{name}.csv')
    ci = fit.conf_int()
    coef = pd.DataFrame({'model': name, 'term': fit.params.index, 'estimate': fit.params.values,
                             'ci_low': ci[0].values, 'ci_high': ci[1].values, 'p': fit.pvalues.values})
    coef['p_fdr'] = np.nan
    grp = coef['term'].str.startswith('C(group')
    if grp.any(): # Benjamini-Hochberg FDR over the group-vs-control terms
        coef.loc[grp, 'p_fdr'] = multipletests(coef.loc[grp, 'p'], method='fdr_bh')[1]
    coef_all.append(coef)
    print(f'\n{name}: n={int(fit.nobs)}, R2={fit.rsquared:.2f}, adj R2={fit.rsquared_adj:.2f}')
    print(anova.round(4).to_string())
pd.concat(coef_all).to_csv(outdir + 'coefficients.csv', index=False)

# all pairwise group comparisons for the volumes, unadjusted for covariates
pairs = []
for y in cols:
    d = df[['group', y]].dropna()
    t = pairwise_tukeyhsd(d[y], d['group'])
    tab = pd.DataFrame(t.summary().data[1:], columns=t.summary().data[0]); tab.insert(0, 'outcome', y); pairs.append(tab)
pd.concat(pairs).to_csv(outdir + 'pairwise_tukey.csv', index=False)
print('\nwritten to', outdir)
