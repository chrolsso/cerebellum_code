import glob
import os
import nibabel as nib



def series_number(path):
    #So that _ph and no _ph gives the same output. 
    name = os.path.basename(path)[:-len('.nii')]
    if name.endswith('_ph'):
        name = name[:-len('_ph')]
    return int(name.split('_')[-1])


def pick_file(subject, pattern):
    #Multiple files match pattern. Let's try taking the most recent one (higher series number). Not sure if correct at the moment. 
    matches = glob.glob(nifti_folder + subject + '/**/' + pattern, recursive=True)
    assert matches, f'{subject}: no file matches {pattern}'
    matches.sort(key=series_number)
    chosen = matches[-1]
    choices[subject, pattern] = [subject, pattern, os.path.basename(chosen), len(matches)]
    return chosen

wdir = '/safe/data/cerebellar-cortex'
nifti_folder = wdir + '/Niftis/'
output_folder = wdir + '/Derivatives/input/'

choices = {}   # which file was used for each subject and pattern, written to selected_files.txt

subjects = sorted(name for name in os.listdir(nifti_folder) if os.path.isdir(nifti_folder + name))

for subject in subjects:
    # the UNI image of the whole-brain scan has the series number of that scan + 4
    whole_brain_series = series_number(pick_file(subject, '*WIP_WB_MP2RAGE*[0-9].nii'))

    # name the pipeline expects, pattern of the raw file, volume number if the file is 4D)
    wanted = [('_anat_inv1',     '*WIP_Cereb_MP2RAGE*[0-9].nii', 0),
              ('_anat_inv2',     '*WIP_Cereb_MP2RAGE*[0-9].nii', 1),
              ('_anat_inv1_ph',  '*WIP_Cereb_MP2RAGE*_ph.nii',   0),
              ('_anat_inv2_ph',  '*WIP_Cereb_MP2RAGE*_ph.nii',   1),
              ('_anat_HEADT1w',  f'*MP2RAGE_UNI*_{whole_brain_series + 4}.nii', None)]

    os.makedirs(f'{output_folder}{subject}/anat', exist_ok=True)

    for suffix, pattern, volume in wanted:
        source = pick_file(subject, pattern)
        target = f'{output_folder}{subject}/anat/{subject}{suffix}.nii.gz'
        if os.path.exists(target):
            continue

        image = nib.load(source)
        if volume is not None:
            data = image.get_fdata(dtype='float32')
            image = nib.Nifti1Image(data[..., volume], image.affine)
        nib.save(image, target)

with open(output_folder + 'selected_files.txt', 'w') as log_file:
    for subject, pattern, used_file, number_of_candidates in choices.values():
        log_file.write(f'{subject}: {used_file}   ({number_of_candidates} candidate(s) for {pattern})\n')

with_several = {row[0] for row in choices.values() if row[3] > 1}
print(len(with_several), 'subject(s) had several candidates, see', output_folder + 'selected_files.txt')
