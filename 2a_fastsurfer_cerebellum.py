import glob
import os
import nibabel
import numpy
from nibabel.processing import resample_from_to
from scipy import ndimage

inpath = '/.../input' # as in the other steps
fastsurfer = '/path/to/the/fastsurfer/output/' # expects <subject>/mri/aseg.auto_noCCseg.mgz
dilate_voxels = 3 # about 2.4 mm on the 0.8 mm whole-head grid, so the dura lies just outside

for anat in sorted(glob.glob(inpath + '*/anat/')):
    subject = anat.replace(inpath, '').split('/')[0]
    aseg_file = fastsurfer + subject + '/mri/aseg.auto_noCCseg.mgz'
    if not os.path.exists(aseg_file):
        print(subject, 'no FastSurfer aseg yet, skipping')
        continue
    aseg = nibabel.load(aseg_file)
    cerebellum = numpy.isin(numpy.asarray(aseg.dataobj), [7, 8, 46, 47]).astype('float32') #cerebellar white matter and cortex, both sides
    head = nibabel.load(anat + subject + '_anat_HEADT1w.nii.gz') # the whole-head grid that step 2 registers the slab to
    on_head = resample_from_to(nibabel.Nifti1Image(cerebellum, aseg.affine), (head.shape, head.affine), order=1).get_fdata() > 0.5
    volume = on_head.sum() * numpy.prod(head.header.get_zooms()[:3]) / 1000
    print(subject, 'cerebellum', round(volume), 'cm3')
    assert 60 < volume < 250, f'{subject}: check the FastSurfer cerebellum, it is {volume:.0f} cm3'
    soft = ndimage.gaussian_filter(ndimage.binary_dilation(on_head, iterations=dilate_voxels).astype('float32'), 1)# soft edge
    out = inpath + subject + '/nighres/coreg/'
    os.makedirs(out, exist_ok=True)
    nibabel.save(nibabel.Nifti1Image(soft, head.affine), out + subject + '_wb_cerebellum_prob.nii.gz')

    # lobular atlas: CerebNet labels, extended to the whole aseg cerebellum by nearest label
    cerebnet_file = fastsurfer + subject + '/mri/cerebellum.CerebNet.nii.gz' # check this name in your FastSurfer output
    if os.path.exists(cerebnet_file):
        labels = numpy.asarray(nibabel.load(cerebnet_file).dataobj).astype('int16')
        aseg_cerebellum = numpy.isin(numpy.asarray(aseg.dataobj), [7, 8, 46, 47])
        assert labels.shape == aseg_cerebellum.shape, f'{subject}: CerebNet and aseg are on different grids'
        nearest = ndimage.distance_transform_edt(labels == 0, return_distances=False, return_indices=True)
        atlas = numpy.where(aseg_cerebellum, labels[tuple(nearest)], 0).astype('float32')
        atlas_on_head = resample_from_to(nibabel.Nifti1Image(atlas, aseg.affine), (head.shape, head.affine), order=0) # order=0 keeps the labels
        nibabel.save(nibabel.Nifti1Image(numpy.asarray(atlas_on_head.dataobj), head.affine), out + subject + '_wb_atlas.nii.gz')
    else:
        print(subject, 'no CerebNet labels, no atlas written')
