#!/usr/bin/env python3
import json, os, shutil, subprocess, sys
import numpy as np, nibabel as nib
from scipy import ndimage as ndi

# ---------------- edit these ----------------
wdir    =  #working directory
T1   = wdir + "BIDS/sub-01/anat/sub-01_T1w.nii.gz" # whole-head T1w (UNI), with skull
SEG  = wdir + "synth/sub-01_synthseg.nii.gz"# SynthSeg labels of that T1w 
TPL  = "/home/fmriprep/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_"
THREADS = "4"
# --------------------------------------------
TASK = nameoftask
BOLD = name of bold file
JSON = corresponding json file 
anatomydir = wdir + "own_pipeline/anat/"     # shared anatomy, computed once
outputdir = wdir + "own_pipeline/" + TASK + "/" # one folder per run
os.makedirs(anatomydir, exist_ok=True)
os.makedirs(outputdir, exist_ok=True)
os.environ["FSLOUTPUTTYPE"] = "NIFTI_GZ"
os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = THREADS

need = ["antsRegistration", "antsApplyTransforms", "N4BiasFieldCorrection",
        "ResampleImageBySpacing", "antsRegistrationSyNQuick.sh", "mcflirt"]
missing = [t for t in need if not shutil.which(t)]
if missing:
    sys.exit("missing tools in this pod: " + " ".join(missing))


def run(cmd):
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def save(arr, ref, path):
    nib.save(nib.Nifti1Image(np.asarray(arr), ref.affine), path)


def have(path):
    return os.path.exists(path)


# ---- 1) anatomy: bias correction, SynthSeg brain mask, T1w -> MNI ---------------------
T1N4 = anatomydir + "T1_n4.nii.gz"
if not have(T1N4):
    run(["N4BiasFieldCorrection", "-d", "3", "-i", T1, "-o", T1N4])
SEGT1 = anatomydir + "seg_in_T1.nii.gz"
if not have(SEGT1):
    run(["antsApplyTransforms", "-d", "3", "-i", SEG, "-r", T1N4, "-t", "identity",
         "-n", "NearestNeighbor", "-o", SEGT1])
T1M, T1B = anatomydir + "T1_mask.nii.gz", anatomydir + "T1_brain.nii.gz"
if not have(T1B):
    t1 = nib.load(T1N4)
    m = ndi.binary_fill_holes(np.rint(nib.load(SEGT1).get_fdata()) > 0)   # labels incl. the thin CSF rim
    save(m.astype("uint8"), t1, T1M)
    save(t1.get_fdata() * m, t1, T1B)
    print("T1 mask covers", round(100 * float(m.mean()), 1), "% of the volume")

MNIB, MNI2 = anatomydir + "mni_brain.nii.gz", anatomydir + "mni_2mm.nii.gz"
if not have(MNI2):
    f1, f2 = TPL + "T1w.nii.gz", TPL + "desc-brain_mask.nii.gz"
    if not (have(f1) and have(f2)) or min(os.path.getsize(f1), os.path.getsize(f2)) < 10000:
        sys.exit("template file missing or an empty placeholder: " + f1 + " / " + f2)
    tp = nib.load(f1)
    save(tp.get_fdata() * (nib.load(f2).get_fdata() > 0), tp, MNIB)
    run(["ResampleImageBySpacing", "3", MNIB, MNI2, "2", "2", "2"])

P1 = anatomydir + "t1_to_mni_"                       # use antsRegistrationSyN.sh instead of ...Quick for final results
if not have(P1 + "1Warp.nii.gz"):
    run(["antsRegistrationSyNQuick.sh", "-d", "3", "-f", MNIB, "-m", T1B, "-o", P1, "-n", THREADS])

# ---- 2) functional: motion correction, reference, EPI mask, EPI -> T1w --------------------
MC = outputdir + "mc.nii.gz"
if not have(MC):
    run(["mcflirt", "-in", BOLD, "-out", outputdir + "mc", "-plots", "-spline_final"])
par = np.loadtxt(outputdir + "mc.par")              # 3 rotations (rad), then 3 translations (mm)
img = nib.load(MC)
data = img.get_fdata(dtype=np.float32)
n = data.shape[-1]
MEAN, MEANN4 = outputdir + "mc_mean.nii.gz", outputdir + "mc_mean_n4.nii.gz"
if not have(MEANN4):
    save(data.mean(-1), img, MEAN)
    run(["N4BiasFieldCorrection", "-d", "3", "-i", MEAN, "-o", MEANN4])

EM = outputdir + "epi_mask.nii.gz"
if not have(EM):
    tmp = outputdir + "t1mask_in_epi.nii.gz"
    run(["antsApplyTransforms", "-d", "3", "-i", T1M, "-r", MEANN4, "-t", "identity",
         "-n", "NearestNeighbor", "-o", tmp])   # header alignment
    ref = nib.load(MEANN4)
    e = ref.get_fdata()
    m = nib.load(tmp).get_fdata() > 0.5
    sig = ndi.gaussian_filter(e, 1.0) > 0.15 * np.percentile(e[m], 98)
    keep = ndi.binary_fill_holes(ndi.binary_closing(m & sig, iterations=2)) & m
    save(keep.astype("uint8"), ref, EM)
    print("EPI mask covers", round(100 * float(keep.mean()), 1), "% of the slab volume")
keep = nib.load(EM).get_fdata() > 0

P2 = outputdir + "own_"
AFF = P2 + "0GenericAffine.mat"
if not have(AFF):
    run(["antsRegistration", "-d", "3", "--float", "1", "-o", "[" + P2 + "," + outputdir + "epi_in_T1w.nii.gz]",
         "-n", "Linear", "-w", "[0.005,0.995]", "-t", "Rigid[0.1]",
         "-m", "MI[" + T1B + "," + MEANN4 + ",1,32,Regular,0.25]",
         "-c", "[1000x500x250x100,1e-7,10]", "-f", "8x4x2x1", "-s", "3x2x1x0vox",
         "-x", "[" + T1M + "," + EM + "]"])

# ---- 3) one ANTs resampling of the whole series into MNI (the last -t is applied first) ------
chain = ["-t", P1 + "1Warp.nii.gz", "-t", P1 + "0GenericAffine.mat", "-t", AFF]
if not have(outputdir + "bold_MNI.nii.gz"):
    run(["antsApplyTransforms", "-d", "3", "-e", "3", "--float", "1", "-i", MC, "-r", MNI2]
        + chain + ["-o", outputdir + "bold_MNI.nii.gz"])
if not have(outputdir + "mean_epi_in_MNI.nii.gz"):
    run(["antsApplyTransforms", "-d", "3", "--float", "1", "-i", MEANN4, "-r", MNI2]
        + chain + ["-o", outputdir + "mean_epi_in_MNI.nii.gz"])

# ---- 4) confounds: motion, FD, DVARS, global/WM/CSF signal, 5 aCompCor components ------------
SEGE = outputdir + "seg_in_epi.nii.gz"
if not have(SEGE):
    run(["antsApplyTransforms", "-d", "3", "-i", SEGT1, "-r", MEANN4, "-t", "[" + AFF + ",1]",
         "-n", "NearestNeighbor", "-o", SEGE])
le = np.rint(nib.load(SEGE).get_fdata()).astype(int)
wm = ndi.binary_erosion(np.isin(le, (2, 41, 7, 46)) & keep)
csf = np.isin(le, (24, 4, 43, 14, 15)) & keep

tt = np.arange(n) - (n - 1) / 2
X = np.vstack([np.ones(n), tt, tt ** 2]).T


def detrend(a):                              # voxels x time
    return (a.T - X @ np.linalg.lstsq(X, a.T, rcond=None)[0]).T


d = np.abs(np.diff(par, axis=0))
fd = np.r_[np.nan, d[:, 3:].sum(1) + 50 * d[:, :3].sum(1)]
ts = data[keep]
dv = np.r_[np.nan, np.sqrt((np.diff(ts, axis=1) ** 2).mean(0))]
roi = detrend(data[wm | csf].astype(np.float64))
roi = roi / np.maximum(roi.std(1, keepdims=True), 1e-6)
U = np.linalg.svd(roi.T, full_matrices=False)[0]
cols = {"trans_x": par[:, 3], "trans_y": par[:, 4], "trans_z": par[:, 5],
        "rot_x": par[:, 0], "rot_y": par[:, 1], "rot_z": par[:, 2],
        "framewise_displacement": fd, "dvars": dv, "global_signal": ts.mean(0),
        "csf": data[csf].mean(0), "white_matter": data[wm].mean(0)}
for i in range(5):
    cols["a_comp_cor_%02d" % i] = U[:, i]
with open(outputdir + "confounds.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for i in range(n):
        f.write("\t".join("n/a" if np.isnan(v[i]) else "%.6g" % v[i] for v in cols.values()) + "\n")
print("FD mean", round(float(np.nanmean(fd)), 3), " max", round(float(np.nanmax(fd)), 3),
      " volumes >0.2 mm:", int((fd > 0.2).sum()), " >0.5 mm:", int((fd > 0.5).sum()), " of", n)

# ---- 5) sensitivity in cerebellar gray vs white matter (on the motion-corrected series) -------
TR = json.load(open(JSON))["RepetitionTime"]
freqs = np.fft.rfftfreq(n, TR)
low = (freqs >= 0.01) & (freqs <= 0.1)
means = {}
for name, ids in (("cerebellar GM", (8, 47)), ("cerebellar WM", (7, 46))):
    sel = np.isin(le, ids) & keep
    if sel.sum() < 50:
        print(name, ": too few voxels inside the mask")
        continue
    t = data[sel].astype(np.float64)
    mean = t.mean(1)
    res = detrend(t)
    sd = res.std(1)
    amp = np.abs(np.fft.rfft(res, axis=1))
    falff = amp[:, low].sum(1) / np.maximum(amp[:, 1:].sum(1), 1e-6)
    means[name] = float(mean.mean())
    print(name, "| voxels", int(sel.sum()), "| mean signal", round(means[name], 1),
          "| detrended tSNR (median)", round(float(np.median(mean / np.maximum(sd, 1e-6))), 1),
          "| temporal SD % of mean (median)", round(float(np.median(100 * sd / mean)), 2),
          "| 0.01-0.1 Hz fraction (median)", round(float(np.median(falff)), 3))
if len(means) == 2:
    g, w = means["cerebellar GM"], means["cerebellar WM"]
    print("tissue contrast in the mean EPI (GM-WM)/mean:", round(100 * (g - w) / ((g + w) / 2), 1), "%")
