import glob
import os
import nibabel
import numpy
from scipy import ndimage

inpath = '/input/'       # as in the other steps
sigma_mm = 1.0        # width of the smoothing of the label boundaries

for coreg in sorted(glb.glob(inpath + '*/nighres/coreg/')):
    subject = coreg.replace(inpath, '').split('/')[0]
    name = lambda part: coreg + subject + '_anat_inv1_s2c10_lcpca-den_qt1map-t1_cr_' + part
    if not (os.path.exists(name('atlas_def-img.nii.gz')) and os.path.exists(name('cbmask_def-img.nii.gz'))):
        print(subject, 'no cropped atlas or mask yet, skipping')
        continue
    atlas_img = nibabel.load(name('atlas_def-img.nii.gz'))
    labels = numpy.rint(numpy.asarray(atlas_img.dataobj)).astype('int16')
    mask = numpy.asarray(nibabel.load(name('cbmask_def-img.nii.gz')).dataobj) > 0.5
    sigma = sigma_mm / float(numpy.mean(atlas_img.header.get_zooms()[:3]))

    # work only inside the box around the mask, which keeps the time and memory down
    pad = int(4 * sigma) + 2
    box = tuple(slice(max(a.min() - pad, 0), a.max() + pad + 1) for a in numpy.nonzero(mask | (labels > 0)))
    lab, msk = labels[box], mask[box]

    # 1. every voxel of the mask gets the label of the nearest atlas voxel
    nearest = ndimage.distance_transform_edt(lab == 0, return_distances=False, return_indices=True)
    filled = numpy.where(msk, lab[tuple(nearest)], 0).astype('int16')

    # 2. smooth the label boundaries: each voxel takes the label with the highest smoothed indicator
    best = numpy.zeros(filled.shape, dtype='float32')
    smooth = numpy.zeros(filled.shape, dtype='int16')
    for value in numpy.unique(filled[filled > 0]):
        g = ndimage.gaussian_filter((filled == value).astype('float32'), sigma)
        better = g > best
        best[better] = g[better]
        smooth[better] = value
    smooth[~msk] = 0

    out = numpy.zeros(labels.shape, dtype='float32')
    out[box] = smooth
    print(subject, 'labelled voxels', int((labels > 0).sum()), '->', int((out > 0).sum()), '| mask voxels', int(mask.sum()))
    nibabel.save(nibabel.Nifti1Image(out, atlas_img.affine, atlas_img.header), name('atlas_smooth_def-img.nii.gz'))
