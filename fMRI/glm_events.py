#!/usr/bin/env python3

import json, os, sys
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.glm.first_level import FirstLevelModel
from nilearn.glm import threshold_stats_img

# ---------------- edit these ----------------
wdir   = ""
TASK   = "run-01" #own pipeline
LOGS   = wdir + ""
EVENTS = LOGS + ".tsv"
PULSES = LOGS + "-.tsv" # "" to skip the check
JSON   = wdir + "Niftis/.json"

MODE   = "task"      #"task"  = every event in ONE regressor vs implicit baseline
                     #(the motor positive control; start here)
                     #or a contrast string over trial_type, e.g. "pattern - random"
T0     = 0.0     #seconds to subtract from the onsets (0 if they start at the first volume)
DUR    = None      #None = use the duration column; a number overrides every duration (s)
FWHM   = 4.0   #smoothing, mm
SHIFTS = [0.0, 40.0]# first = real, the rest are negative controls
# --------------------------------------------
outputdir = wdir + "own_pipeline/sub-01/" + TASK + "/"
BOLD=outputdir + "mc.nii.gz"
MASK=outputdir + "epi_mask.nii.gz"
CONF=outputdir + "confounds.tsv"
LAB = outputdir + "seg_in_epi.nii.gz"
for p in (BOLD, MASK, CONF, EVENTS, JSON):
    if not os.path.exists(p):
        sys.exit("missing: " + p)

TR = json.load(open(JSON))["RepetitionTime"]
n = nib.load(BOLD).shape[-1]
LEN = n * TR
print("volumes %d | TR %.3f s | scan length %.1f s" % (n, TR, LEN))

if PULSES and os.path.exists(PULSES):
    npulse = sum(1 for _ in open(PULSES)) - 1          # assumes one header row
    print("pulses in the log: %d (volumes: %d)" % (npulse, n))
    if abs(npulse - n) > 2:
        print("   NOTE: the log does not cover the whole scan; check which end is missing")

# ---- events ----------------------------------------------------------------
ev = pd.read_csv(EVENTS, sep="\t")
print("columns:", list(ev.columns))
if "onset" not in ev.columns:
    sys.exit("no 'onset' column, check the events file")
if "duration" not in ev.columns:
    ev["duration"] = 0.0
if "trial_type" not in ev.columns:
    ev["trial_type"] = "task"
ev = ev.dropna(subset=["onset"]).copy()
ev["trial_type"] = ev["trial_type"].astype(str)
print(ev.head(3).to_string())
print("rows:", len(ev), " trial types:", dict(ev["trial_type"].value_counts()))

if ev["onset"].max() > 5 * LEN:  # onsets look like milliseconds
    print("onsets exceed 5x the scan length -> treating them as milliseconds")
    ev["onset"] = ev["onset"] / 1000.0
    ev["duration"] = ev["duration"] / 1000.0
ev["onset"] = ev["onset"] - T0
if DUR is not None:
    ev["duration"] = float(DUR)
print("first onset %.1f s | last onset %.1f s | scan ends %.1f s" % (ev["onset"].min(), ev["onset"].max(), LEN))
drop = (ev["onset"] < 0) | (ev["onset"] >= LEN)
if drop.any():
    print("   dropping %d events outside the scan" % int(drop.sum()))
    ev = ev[~drop]
if len(ev) == 0:
    sys.exit("no events fall inside the scan: check T0 and the onset units")

cover = float(np.clip(ev["duration"], 0.5, None).sum()) / LEN
print("events cover roughly %.0f%% of the scan" % (100 * cover))
if cover > 0.9 and MODE == "task":
    print(" WARNING: almost no implicit baseline, so 'task vs rest' will be weak.")
    print("If the task ran in blocks, set MODE to a contrast between conditions instead.")

if MODE == "task":
    ev = ev.assign(trial_type="task")
    contrast = "task"
else:
    contrast = MODE
ev = ev[["onset", "duration", "trial_type"]]

# ---- confounds -------------------------------------------------------------
conf = pd.read_csv(CONF, sep="\t")
use = [c for c in ("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z") if c in conf.columns]
use += [c for c in conf.columns if c.startswith("a_comp_cor")]
conf = conf[use].fillna(0.0)
if len(conf) != n:
    sys.exit("confounds has %d rows but the BOLD has %d volumes" % (len(conf), n))
print("confound regressors:", len(use))

# ---- fit, real timing then the negative controls ----------------------------
inmask = nib.load(MASK).get_fdata() > 0
lab = np.rint(nib.load(LAB).get_fdata()).astype(int) if os.path.exists(LAB) else None
if lab is None:
    print("no", LAB, "-> per-tissue counts skipped")

res = {}
for sh in SHIFTS:
    e = ev.copy()
    e["onset"] = e["onset"] + sh
    e = e[(e["onset"] >= 0) & (e["onset"] < LEN)]
    model = FirstLevelModel(t_r=TR, hrf_model="spm", drift_model="cosine", high_pass=1 / 128.0,
                            noise_model="ar1", smoothing_fwhm=FWHM, mask_img=MASK, minimize_memory=False)
    model.fit(BOLD, events=e, confounds=conf)
    z = model.compute_contrast(contrast, output_type="z_score")
    nib.save(z, outputdir + "glm_z_shift%d.nii.gz" % int(sh))
    zd = np.nan_to_num(z.get_fdata())
    hit = (zd > 3.1) & inmask
    r = {"max_z": float(zd[inmask].max()), "hit": int(hit.sum())}
    line = "shift %4.0f s | max z %5.2f | voxels z>3.1: %6d of %d" % (sh, r["max_z"], r["hit"], int(inmask.sum()))
    if lab is not None:
        r["gm"] = int((hit & np.isin(lab, (8, 47))).sum())
        r["wm"] = int((hit & np.isin(lab, (7, 46))).sum())
        line += " | cerebellar GM %d, WM %d" % (r["gm"], r["wm"])
    res[sh] = r
    print(line, flush=True)

# ---- verdict ----------------------------------------------------------------
real = res[SHIFTS[0]]
ctrl = [res[s] for s in SHIFTS[1:]]
if ctrl:
    worst = max(c["hit"] for c in ctrl)
    print("\nreal timing: %d suprathreshold voxels | best negative control: %d" % (real["hit"], worst))
    if real["hit"] > 3 * max(worst, 1):
        print("=> the run carries a clear task response")
    elif real["hit"] > worst:
        print("=> weak: a response may be there, but not convincingly above the control")
    else:
        print("=> no detectable response with this contrast; check the timing before blaming the data")

try:
    thr = threshold_stats_img(nib.load(outputdir + "glm_z_shift0.nii.gz"), mask_img=MASK,
                              alpha=0.05, height_control="fdr", two_sided=False)[1]
    print("FDR q<0.05 threshold (z): %.2f" % float(thr))
except Exception as exc:
    print("FDR threshold not computed:", exc)
print("z-maps written to", outputdir)
