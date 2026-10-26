import glob
import os
import re
import nibabel as nib
import numpy as np
import pandas as pd

#initial settings
excel_file = '/path/to/the/clinical.xlsx'
sheet = 0
col_id = 'ID'  # First column name
col_date = 'DATUM MR' #scan date, e.g. 2025-08-12
col_initials = 'INITIALER' # initials column
col_group = 'DIAGNOS'  # diagnosis column
disease_values = []  # column entries for diagnosis every other entry counts as Control
keep_columns = []  # variable names to keep. 
fastsurfer = '/fastsurfer/output/'# <subject>/mri/aseg.auto_noCCseg.mgz
out_file = 'output/clinical_merged.csv'
key_file = 'output/subject_key.csv'  
#%%

clean = lambda name: re.sub(r'\W+', '_', str(name)).strip('_').lower() # column names usable in formulas

table = pd.read_excel(excel_file, sheet_name=sheet)
table['_date'] = pd.to_datetime(table[col_date], errors='coerce').dt.strftime('%Y%m%d')
table['_ini'] = table[col_initials].astype(str).str.replace(r'[^A-Za-z]', '', regex=True).str.upper()

#volumes from the FastSurfer aseg (cm3)
#matches cerebellum labels to the real measurements
labels = {'cerebellum_wm': [7, 46], 'cerebellum_cortex': [8, 47], 'cerebellum_left': [7, 8], 'cerebellum_right': [46, 47]}
rows, problems = [], []

for aseg_file in sorted(glob.glob(fastsurfer + '*/mri/aseg.auto_noCCseg.mgz')):


    subject = aseg_file.replace(fastsurfer, '').split('/')[0]

    m = re.match(r'^(\d{8})_([A-Za-z]+)', subject)  # split the subject name into the 8-digit date and the initials
    if not m:
        problems.append(f'{subject}: name is not DATE_INITIALS')
        continue
    
    date, ini = m.group(1), m.group(2).upper()
    hit = table[(table['_date'] == date) & (table['_ini'] == ini)] #

    if len(hit) != 1:
        problems.append(f'{subject}: {len(hit)} matching rows in the Excel file')
        continue

    seg = nib.load(aseg_file)
    data = np.asarray(seg.dataobj)
    voxel = float(np.prod(seg.header.get_zooms()[:3])) / 1000
    vol = {k: np.isin(data, v).sum() * voxel for k, v in labels.items()}
    vol['cerebellum_total'] = vol['cerebellum_wm'] + vol['cerebellum_cortex']
    vol['brain_total'] = (data > 0).sum() * voxel # all labelled voxels of the aseg
    row = hit.iloc[0]
    group = str(row[col_group]).strip() if str(row[col_group]).strip() in disease_values else 'Control'
    rows.append({'subject': subject, 'id': str(row[col_id]).replace(' ', '_'), 'group': group,
                 **{clean(c): row[c] for c in keep_columns}, **{k: round(v, 2) for k, v in vol.items()}})

result = pd.DataFrame(rows)
if len(result):
    result.drop(columns='subject').to_csv(out_file, index=False) # no subject names in the file that goes to plotting
    result[['subject', 'id']].to_csv(key_file, index=False)
    print(result.groupby('group').size().rename('n'))
print(len(rows), 'subjects merged,', len(problems), 'problems')
for p in problems: print('  ', p)
