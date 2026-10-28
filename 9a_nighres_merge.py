import glob
import numpy as np
import pandas as pd
import nibabel as nib

#initial settings
root = 'input/'  
outdir = 'output/'
key_file = outdir + 'subject_key.csv' #Take these from 8a
merged = outdir + 'clinical_merged.csv' 
qc_file = root + 'qc_summary.csv'#From 6_qc.py
age_col, sex_col = 'ålder', 'kön'#Names of these columns in clinical_merged.csv
exclude = []#IDs to leave out, e.g. the subject with motion artefacts
thickness_unit_mm = None #Double check this what it should be
min_voxels = 200 #Regions with fewer cortex voxels get no thickness
keep_qc = ['cb_mask_cm3', 'cb_mask_components', 'wm_rawprior_components', 'wm_object_components',
           'share_class1', 'share_class2', 'share_class3', 'share_class4', 'cortex_cm3', 'thickness_median_mm', 'thickness_p95_mm']
#%%

def load(subject, folder, pattern):
    files = glob.glob(f'{root}{subject}/nighres/{folder}/{pattern}')
    assert len(files) == 1, (subject, pattern, files)
    return nib.load(files[0])

key = pd.read_csv(key_file)
info = pd.read_csv(merged).rename(columns={age_col: 'age', sex_col: 'sex'})[['id', 'group', 'status', 'age', 'sex', 'tiv']]

rows = []
for subject, id_ in key.itertuples(index=False):
    if id_ in exclude: continue
    try:
        thick_img = load(subject, 'surfaces', '*cruise-thick.nii.gz')
        cortex = load(subject, 'surfaces', '*cruise-cortex.nii.gz').get_fdata() > 0.5
        gm = load(subject, 'tissues', '*_msk_gmrawprior.nii.gz').get_fdata()
        wm = load(subject, 'tissues', '*_msk_wmrawprior.nii.gz').get_fdata()
        thick = thick_img.get_fdata()
        zooms = thick_img.header.get_zooms()[:3]
        scale = thickness_unit_mm or float(np.mean(zooms))
        voxel_cm3 = float(np.prod(zooms)) / 1000
        regions = {'all': load(subject, 'coreg', '*_cr_cbmask_def-img.nii.gz').get_fdata() > 0.5,
                   'left': load(subject, 'coreg', '*_cr_cbmask_l_def-img.nii.gz').get_fdata() > 0.5,
                   'right': load(subject, 'coreg', '*_cr_cbmask_r_def-img.nii.gz').get_fdata() > 0.5}
    except Exception as err:
        print('Skipping', subject, ':', err)
        continue
    try:
        atlas = np.rint(load(subject, 'coreg', '*_cr_atlas_smooth_def-img.nii.gz').get_fdata()).astype(int)
        for lab in np.unique(atlas[atlas > 0]): regions[f'lobule_{lab}'] = atlas == lab
    except Exception as err:
        print('No lobules for', subject, ':', err)

    for region, mask in regions.items():
        t = thick[cortex & mask & (thick > 0)] * scale
        ok = t.size >= min_voxels
        values = {'thickness_median_mm': np.median(t) if ok else np.nan, 'thickness_mean_mm': t.mean() if ok else np.nan,
                  'thickness_p10_mm': np.percentile(t, 10) if ok else np.nan, 'thickness_p90_mm': np.percentile(t, 90) if ok else np.nan,
                  'cortex_volume_cm3': (cortex & mask).sum() * voxel_cm3, 'gm_volume_cm3': (gm * mask).sum() * voxel_cm3,
                  'wm_volume_cm3': (wm * mask).sum() * voxel_cm3, 'region_volume_cm3': mask.sum() * voxel_cm3}
        for measure, value in values.items():
            rows.append({'id': id_, 'region': region, 'measure': measure, 'value': float(value)})
    print(subject, len(regions), 'regions')

assert rows, 'No subject had all the files'
long = pd.DataFrame(rows).merge(info, on='id', how='left')       # ids only, no subject names in the output
long.to_csv(outdir + 'nighres_long.csv', index=False)

qc = pd.read_csv(qc_file).merge(key, on='subject')
qc[['id'] + [c for c in keep_qc if c in qc.columns]].to_csv(outdir + 'nighres_qc.csv', index=False)
print(long['id'].nunique(), 'subjects,', long['region'].nunique(), 'regions written to', outdir)
