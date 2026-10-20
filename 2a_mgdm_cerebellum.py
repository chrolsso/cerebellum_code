import glob
import nibabel
import numpy
import nighres
from scipy import ndimage

inpath = '/data/safe/cerebellar-cortex-CO2026/Derivatives/input/'

for work in glob.glob(inpath + '*/tiv/'):
    subject = work.replace(inpath, '').split('/')[0]
    strip = nighres.brain.mp2rage_skullstripping(second_inversion=work + 'inv2.nii.gz', t1_map=work + 't1map.nii.gz',
                                                 save_data=True, output_dir=work)
    brain = nibabel.load(strip['brain_mask']).get_fdata() > 0.5
    t1_img = nibabel.load(strip['t1map_masked'])
    t1 = t1_img.get_fdata(dtype='float32')

    hole = brain & (t1 <= 0)
    if hole.sum() > 0.005 * brain.sum():   # the scanner T1 map has zeros inside the brain: fill them from the UNI
        uni = nibabel.as_closest_canonical(nibabel.load(inpath + subject + '/anat/' + subject + '_anat_HEADT1w.nii.gz')).get_fdata(dtype='float32')
        ok = brain & (t1 > 0)
        edges = numpy.unique(numpy.percentile(uni[ok], numpy.linspace(0, 100, 201)))
        bins = numpy.clip(numpy.digitize(uni[ok], edges) - 1, 0, len(edges) - 2) + 1
        median = numpy.array(ndimage.median(t1[ok], bins, numpy.arange(1, len(edges))))
        centres = (edges[:-1] + edges[1:]) / 2
        valid = numpy.isfinite(median)
        t1[hole] = numpy.interp(uni[hole], centres[valid], median[valid])
    nibabel.save(nibabel.Nifti1Image(t1, t1_img.affine), work + 't1map_filled.nii.gz')

    mgdm = nighres.brain.mgdm_segmentation(contrast_image1=work + 't1map_filled.nii.gz', contrast_type1='T1map7T',
                                           save_data=True, output_dir=work, file_name=subject + '_wb')
    cereb = nighres.brain.extract_brain_region(segmentation=mgdm['segmentation'], levelset_boundary=mgdm['distance'],
                                               maximum_membership=mgdm['memberships'], maximum_label=mgdm['labels'],
                                               extracted_region='cerebellum', save_data=True, output_dir=work, file_name=subject + '_wb')
    gm, wm = nibabel.load(cereb['region_proba']), nibabel.load(cereb['inside_proba'])
    cerebellum = (gm.get_fdata(dtype='float32') + wm.get_fdata(dtype='float32')) > 0.5
    volume = cerebellum.sum() * numpy.prod(gm.header.get_zooms()[:3]) / 1000
    print(subject, 'cerebellum', round(volume), 'cm3')
    assert 60 < volume < 250, f'{subject}: MGDM cerebellum is {volume:.0f} cm3, check it'
    soft = ndimage.gaussian_filter(ndimage.binary_dilation(cerebellum, iterations=3).astype('float32'), 1)   # about 2 mm margin, soft edge
    nibabel.save(nibabel.Nifti1Image(soft, gm.affine), inpath + subject + '/nighres/coreg/' + subject + '_wb_cerebellum_prob.nii.gz')
