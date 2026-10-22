import glob
import os
import nibabel
import nighres
import numpy
from nibabel.processing import resample_from_to
from scipy import ndimage

inpath = '/Postprocessed/input/'       # as in the other steps; needs <subject>/anat/<subject>_anat_MPRAGE.nii.gz
fastsurfer = '/path/to/the/fastsurfer/output/'                          # expects <subject>/mri/aseg.auto_noCCseg.mgz, run on that MPRAGE
dilate_voxels = 3                                                       # about 2.4 mm on the 0.8 mm whole-head grid, so the dura lies just outside

for anat in sorted(glob.glob(inpath + '*/anat/')):
    subject = anat.replace(inpath, '').split('/')[0]
    aseg_file = fastsurfer + subject + '/mri/aseg.auto_noCCseg.mgz'
    if not os.path.exists(aseg_file):
        print(subject, 'no FastSurfer aseg yet, skipping')
        continue
    work = inpath + subject + '/nighres/coreg/'
    os.makedirs(work, exist_ok=True)
    mprage_file = anat + subject + '_anat_MPRAGE.nii.gz'
    head_file = anat + subject + '_anat_HEADT1w.nii.gz'

    # 1. cerebellum from the aseg, on the MPRAGE's own grid (headers only, no registration)
    aseg = nibabel.load(aseg_file)
    mprage = nibabel.load(mprage_file)
    cerebellum = numpy.isin(numpy.asarray(aseg.dataobj), [7, 8, 46, 47]).astype('float32')
    on_mprage = resample_from_to(nibabel.Nifti1Image(cerebellum, aseg.affine), (mprage.shape, mprage.affine), order=1).get_fdata() > 0.5
    cerebellum_file = work + subject + '_mprage_cerebellum.nii.gz'
    nibabel.save(nibabel.Nifti1Image(on_mprage.astype('float32'), mprage.affine, mprage.header), cerebellum_file)

    # 2. rigid registration of the MPRAGE to the whole-head UNI (the head may have moved between the two scans)
    reg = nighres.registration.embedded_antspy(source_image=mprage_file, target_image=head_file,
                                               run_rigid=True, rigid_iterations=1000, run_affine=False, run_syn=False,
                                               cost_function='MutualInformation', interpolation='Linear',
                                               save_data=True, overwrite=False, output_dir=work, file_name=subject + '_mprage-to-head')

    # 3. the cerebellum onto the HEADT1w grid
    moved = nighres.registration.apply_coordinate_mappings(cerebellum_file, mapping1=reg['mapping'], interpolation='linear', padding='zero',
                                                           save_data=True, overwrite=False, output_dir=work, file_name=subject + '_cerebellum_on_head.nii.gz')['result']
    head = nibabel.load(head_file)
    on_head = nibabel.load(moved).get_fdata() > 0.5
    volume = on_head.sum() * numpy.prod(head.header.get_zooms()[:3]) / 1000
    print(subject, 'cerebellum', round(volume), 'cm3')
    assert 60 < volume < 250, f'{subject}: check the FastSurfer cerebellum, it is {volume:.0f} cm3'
    soft = ndimage.binary_dilation(on_head, iterations=dilate_voxels).astype('float32')     # the linear resampling in step 2 gives the soft edge
    nibabel.save(nibabel.Nifti1Image(soft, head.affine), work + subject + '_wb_cerebellum_prob.nii.gz')

    # 4. lobular atlas: CerebNet labels, extended to the whole aseg cerebellum by nearest label, moved with the same rigid mapping
    cerebnet_file = fastsurfer + subject + '/mri/cerebellum.CerebNet.nii.gz'                  # check this name in your FastSurfer output
    if os.path.exists(cerebnet_file):
        cerebnet = nibabel.load(cerebnet_file)
        print(subject, 'CerebNet grid', cerebnet.shape, 'aseg grid', aseg.shape)
        labels = numpy.rint(numpy.asarray(resample_from_to(cerebnet, (aseg.shape, aseg.affine), order=0).dataobj)).astype('int16')   # onto the aseg grid, by the headers
        nearest = ndimage.distance_transform_edt(labels == 0, return_distances=False, return_indices=True)
        atlas = numpy.where(numpy.isin(numpy.asarray(aseg.dataobj), [7, 8, 46, 47]), labels[tuple(nearest)], 0).astype('float32')
        atlas_on_mprage = resample_from_to(nibabel.Nifti1Image(atlas, aseg.affine), (mprage.shape, mprage.affine), order=0)           # order=0 keeps the labels
        atlas_file = work + subject + '_mprage_atlas.nii.gz'
        nibabel.save(nibabel.Nifti1Image(numpy.asarray(atlas_on_mprage.dataobj).astype('float32'), mprage.affine, mprage.header), atlas_file)
        moved_atlas = nighres.registration.apply_coordinate_mappings(atlas_file, mapping1=reg['mapping'], interpolation='nearest', padding='zero',
                                                                     save_data=True, overwrite=False, output_dir=work, file_name=subject + '_atlas_on_head.nii.gz')['result']
        nibabel.save(nibabel.Nifti1Image(numpy.rint(nibabel.load(moved_atlas).get_fdata()).astype('float32'), head.affine), work + subject + '_wb_atlas.nii.gz')
    else:
        print(subject, 'no CerebNet labels, no atlas written')
