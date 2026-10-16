import glob
import os
import nibabel
import numpy
from scipy import ndimage

inpath = '/data/safe/cerebellar-cortex-CO2026/Derivatives/input/'       # as in the other steps
radius_mm = 10        # the mask may grow at most this far
csf_factor = 1.5      # growth stops at voxels with T1 above this times the median T1 inside the mask

for coreg in glob.glob(inpath + '*/nighres/coreg/'):
    subject = coreg.replace(inpath, '').split('/')[0]
    name = lambda part: coreg + subject + '_anat_inv1_s2c10_lcpca-den_qt1map-t1_cr_' + part
    if not os.path.exists(name('cbmask_def-img.nii.gz')):
        print(subject, 'no cropped mask yet, skipping')
        continue
    t1 = nibabel.load(name('def-img.nii.gz')).get_fdata(dtype='float32')
    soft = nibabel.load(name('cbmask_def-img.nii.gz'))
    core = soft.get_fdata() > 0.5
    voxel_mm = float(numpy.mean(soft.header.get_zooms()[:3]))
    save = lambda data, part: nibabel.save(nibabel.Nifti1Image(data.astype('float32'), soft.affine, soft.header), name(part))

    # grow the mask through tissue only, at most radius_mm from the original mask
    near = ndimage.distance_transform_edt(~core) <= radius_mm / voxel_mm
    tissue = (t1 > 0) & (t1 < csf_factor * numpy.median(t1[core & (t1 > 0)]))
    grown = ndimage.binary_fill_holes(ndimage.binary_propagation(core, mask=(tissue & near) | core))
    print(subject, 'mask grew from', round(core.sum() * voxel_mm ** 3 / 1000), 'to', round(grown.sum() * voxel_mm ** 3 / 1000), 'cm3')
    save(ndimage.gaussian_filter(grown.astype('float32'), 2), 'cbmask_grown_def-img.nii.gz')     # soft edge, like the original

    # left / right: the side of the nearest original voxel
    left = nibabel.load(name('cbmask_l_def-img.nii.gz')).get_fdata() > 0.5
    right = nibabel.load(name('cbmask_r_def-img.nii.gz')).get_fdata() > 0.5
    nearest_side = ndimage.distance_transform_edt(~(left | right), return_distances=False, return_indices=True)
    is_left = left[tuple(nearest_side)]
    save(grown & is_left, 'cbmask_grown_l_def-img.nii.gz')
    save(grown & ~is_left, 'cbmask_grown_r_def-img.nii.gz')

    # lobular atlas: every grown voxel gets the label of its nearest atlas voxel
    if os.path.exists(name('atlas_def-img.nii.gz')):
        labels = numpy.asarray(nibabel.load(name('atlas_def-img.nii.gz')).dataobj)
        nearest_label = ndimage.distance_transform_edt(labels == 0, return_distances=False, return_indices=True)
        save(numpy.where(grown, labels[tuple(nearest_label)], 0), 'atlas_grown_def-img.nii.gz')
