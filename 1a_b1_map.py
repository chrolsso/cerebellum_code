import glob
import os
import shutil
import nibabel as nib
import nighres
import numpy as np
import scipy.ndimage
from nilearn.image import resample_to_img

# standalone: run after 0_select_files.py and before step 1 (denoising and T1 mapping)
wdir = '/safe/data/cerebellar-cortex'
nifti_folder = wdir + '/Niftis/'
inpath = wdir + '/Derivatives/input/'   # as in the other steps

# per subject (name of its folder in Niftis), the series number of the DREAMB1_combined file that goes with the
# cerebellar scan, e.g. {'SUBJECT': 503}. Its DREAM_anat image is assumed to be that series + 1.
b1_series = {}


def find(subject, pattern):
    matches = glob.glob(nifti_folder + subject + '/**/' + pattern, recursive=True)
    assert len(matches) == 1, f'{subject}: {len(matches)} files match {pattern}'
    return matches[0]


for subject, number in b1_series.items():
    anat = inpath + subject + '/anat/'
    inv2 = anat + subject + '_anat_inv2.nii.gz'      # second inversion of the cerebellar scan from step 0, the target space
    b1ratio = anat + subject + '_b1ratio.nii.gz'     # the map that step 1 uses as B1_map
    if os.path.exists(b1ratio):
        continue
    dream = find(subject, f'*DREAMB1_combined*_{number}.nii')    # DREAM B1 map, in percent of the nominal flip angle
    dream_anat = find(subject, f'*DREAM_anat*_{number + 1}.nii')  # anatomical image of the same DREAM scan
    work = inpath + subject + '/nighres/b1/'
    os.makedirs(work, exist_ok=True)

    # percent to ratio, voxels without a value get the nearest valid one
    img = nib.load(dream)
    ratio = np.asarray(img.dataobj, dtype=np.float32) / 100
    nearest = scipy.ndimage.distance_transform_edt(ratio <= 0, return_distances=False, return_indices=True)
    ratio_file = work + subject + '_dream_b1ratio.nii.gz'
    nib.save(nib.Nifti1Image(ratio[tuple(nearest)], img.affine), ratio_file)

    # the DREAM scan was acquired separately, so move its map onto the cerebellar scan with a rigid registration
    reg = nighres.registration.embedded_antspy(dream_anat, inv2, run_rigid=True, rigid_iterations=1000,
                run_affine=False, run_syn=False, cost_function='MutualInformation', interpolation='Linear',
                save_data=True, overwrite=False, output_dir=work, file_name=subject + '_dream_anat.nii.gz')
    registered = nighres.registration.apply_coordinate_mappings(ratio_file, mapping1=reg['mapping'],
                interpolation='linear', padding='closest', save_data=True, overwrite=False,
                output_dir=work, file_name=subject + '_dream_b1ratio_reg.nii.gz')['result']

    # check against the map resampled by the scanner coordinates only, and use that one if the registration went wrong
    by_header = work + subject + '_dream_b1ratio_header.nii.gz'
    if not os.path.exists(by_header):
        resample_to_img(ratio_file, inv2, interpolation='continuous', force_resample=True, copy_header=True).to_filename(by_header)
    a = nib.load(registered).get_fdata()[::4, ::4, ::4].ravel()
    h = nib.load(by_header).get_fdata()[::4, ::4, ::4].ravel()
    r = np.corrcoef(a, h)[0, 1]
    use = registered if r >= 0.8 else by_header
    shutil.copy(use, b1ratio)

    # B1 ratio over the tissue of the slab, expected to be roughly 0.5 to 1.3 in the cerebellum
    brain = nib.load(inv2).get_fdata(dtype=np.float32)
    tissue = nib.load(b1ratio).get_fdata(dtype=np.float32)[brain > 0.3 * np.percentile(brain, 99)]
    print(subject, f'B1 ratio: registered vs scanner-coordinate resampling r = {r:.2f}{"" if r >= 0.8 else " (using scanner coordinates)"}, median {np.median(tissue):.2f}')
