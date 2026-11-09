import glob
import nighres
import nibabel as nib
import numpy as np
import pandas as pd

#initial settings
root = 'input/'               # folder with one sub-folder per subject, as in 9a
subject = 'SUBJECT'           # folder name of the subject to test
subject_id = 'ID'             # id used in the output files
n_layers = 8                  # layers across the cortex, which gives n_layers + 1 sampled surfaces
min_voxels = 200              # regions with fewer cortex voxels are skipped
n_export = 3000               # single-voxel profiles exported per region, for plotting in MATLAB
outdir = 'output/laminar/'
#%%

base = root + subject + '/nighres/'


def find(folder, pattern):
    files = glob.glob(base + folder + '/' + pattern)
    assert len(files) == 1, f'{pattern}: {len(files)} files found'
    return files[0]


gwb = find('surfaces', '*cruise-gwb.nii.gz')          # WM/GM level set from step 5
cgb = find('surfaces', '*cruise-cgb.nii.gz')          # GM/CSF level set from step 5
cortex = find('surfaces', '*cruise-cortex.nii.gz')
t1 = find('coreg', '*qt1map-t1_cr_def-img.nii.gz')    # the T1 map the segmentation was made on
atlas = find('coreg', '*_cr_atlas_smooth_def-img.nii.gz')
print('inner:', gwb, '\nouter:', cgb, '\nintensity:', t1)

# equivolume layers between the two surfaces: depth 0 at the WM/GM border and 1 at the GM/CSF border
layering = nighres.laminar.volumetric_layering(gwb, cgb, n_layers=n_layers, method='volume-preserving', layer_dir='outward',
                                               save_data=True, overwrite=False, output_dir=outdir, file_name=f'{subject_id}_layering')
# T1 sampled on every layer boundary
sampling = nighres.laminar.profile_sampling(layering['boundaries'], t1, save_data=True, overwrite=False,
                                            output_dir=outdir, file_name=f'{subject_id}_profile')

profile = nib.load(sampling['result']).get_fdata()
depth_img = nib.load(layering['depth']).get_fdata()
cortex = nib.load(cortex).get_fdata() > 0.5
labels = np.round(nib.load(atlas).get_fdata()).astype(int)
print('profile shape', profile.shape, '| cortex voxels', cortex.sum(), '| depth range in cortex', np.round(np.percentile(depth_img[cortex], [1, 99]), 2))
assert profile.shape[-1] == n_layers + 1, 'unexpected number of sampled surfaces, check the layering output'

depth = np.linspace(0, 1, profile.shape[-1])
regions = {'all': cortex}
for label in np.unique(labels[cortex]):
    if label > 0 and (cortex & (labels == label)).sum() >= min_voxels:
        regions[f'lobule_{label}'] = cortex & (labels == label)

rows = []
for name, region in regions.items():
    values = profile[region]                            # voxels x layers
    values = values[(values > 0).all(axis=1)]           # voxels with a valid sample on every surface
    mean, sd = values.mean(axis=0), values.std(axis=0)
    gradient = np.gradient(mean, depth)
    for i in range(len(depth)):
        rows.append({'id': subject_id, 'region': name, 'depth': depth[i], 'mean_t1': mean[i], 'sd_t1': sd[i], 'gradient': gradient[i], 'n_voxels': len(values)})
long = pd.DataFrame(rows)
long.to_csv(outdir + f'{subject_id}_laminar_profile.csv', index=False)

# where the steepest change is, away from the two outer surfaces where partial volume dominates
for name in regions:
    sub = long[(long['region'] == name) & (long['depth'] > 0.1) & (long['depth'] < 0.9)]
    peak = sub.loc[sub['gradient'].abs().idxmax()]
    print(f'{name:>12}: n={int(peak["n_voxels"]):>8}, steepest change at depth {peak["depth"]:.2f}, gradient {peak["gradient"]:.2f}')

# a random sample of single-voxel profiles per region (with voxel indices) for plotting in MATLAB
rng = np.random.default_rng(0)
samples = []
for name, region in regions.items():
    ijk = np.argwhere(region & (profile > 0).all(axis=-1))
    ijk = ijk[rng.choice(len(ijk), min(n_export, len(ijk)), replace=False)]
    table = pd.DataFrame(profile[tuple(ijk.T)], columns=[f't1_{i}' for i in range(profile.shape[-1])])
    table.insert(0, 'k', ijk[:, 2]); table.insert(0, 'j', ijk[:, 1]); table.insert(0, 'i', ijk[:, 0])
    table.insert(0, 'region', name); table.insert(0, 'id', subject_id)
    samples.append(table)
pd.concat(samples).to_csv(outdir + f'{subject_id}_laminar_voxels.csv', index=False)
print('written to', outdir)
