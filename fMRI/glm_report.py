#!/usr/bin/env python3
import json, os, sys
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nilearn.glm.first_level import FirstLevelModel
from nilearn.plotting import plot_design_matrix, plot_stat_map
from nilearn.reporting import get_clusters_table

# ---------------- edit these ----------------
wdir   = ""
TASK   = "run-01" #own pipeline
LOGS   = wdir + ""
EVENTS = LOGS + ".tsv"
PULSES = LOGS + "-.tsv" # "" to skip the check
JSON   = wdir + "Niftis/.json"

MODE   = "task"   #"task" = every event in ONE regressor vs implicit baseline or a contrast string over trial_type, e.g. "high - low"
REST   = ("rest", "fixation")  #trial types kept OUT of the task regressor, so they
                               #stay in the baseline. Set to () to model everything.
T0     = 0.0    #seconds to subtract from the onsets (the first scanner pulse)
DUR    = None   #None = use the duration column; a number overrides every duration (s)
FWHM   = 4.0    #smoothing, mm
THR    = 3.1    #z threshold for the pictures and the cluster table
# --------------------------------------------
outputdir = wdir + "own_pipeline/sub-01/" + TASK + "/"
REP  = outputdir + "report/"
BOLD = outputdir + "mc.nii.gz"
MASK = outputdir + "epi_mask.nii.gz"
CONF = outputdir + "confounds.tsv"
LAB  = outputdir + "seg_in_epi.nii.gz"
BG   = outputdir + "mc_mean_n4.nii.gz"
for p in (BOLD, MASK, CONF, EVENTS, JSON):
    if not os.path.exists(p):
        sys.exit("missing: " + p)
os.makedirs(REP, exist_ok=True)

TR = json.load(open(JSON))["RepetitionTime"]
n = nib.load(BOLD).shape[-1]
LEN = n * TR
print("volumes %d | TR %.3f s | scan length %.1f s" % (n, TR, LEN))

# ---- the pulse log tells you what T0 should be -----------------------------
if PULSES and os.path.exists(PULSES):
    pl = pd.read_csv(PULSES, sep="\t")
    col = [c for c in pl.columns if "onset" in c.lower() or "time" in c.lower()]
    if col:
        t = pd.to_numeric(pl[col[0]], errors="coerce").dropna().values
        print("pulses %d | first %.4f s | last %.1f s | spacing %.4f s | last/TR = %.1f volumes"
              % (len(t), t[0], t[-1], np.median(np.diff(t)), t[-1] / TR))
        print("   -> T0 should be the FIRST pulse time: %.4f" % t[0])
        if abs(t[0] - T0) > 0.05:
            print("   WARNING: T0 is set to %.4f, which does not match" % T0)
        if t[-1] > 1.5 * LEN:
            print("   NOTE: the pulses run well past this scan; the log probably spans more than one run")

# ---- events ----------------------------------------------------------------
ev = pd.read_csv(EVENTS, sep="\t")
if "duration" not in ev.columns:
    ev["duration"] = 0.0
if "trial_type" not in ev.columns:
    ev["trial_type"] = "task"
ev = ev.dropna(subset=["onset"]).copy()
ev["trial_type"] = ev["trial_type"].astype(str)
if ev["onset"].max() > 5 * LEN:
    print("onsets exceed 5x the scan length -> treating them as milliseconds")
    ev["onset"] = ev["onset"] / 1000.0
    ev["duration"] = ev["duration"] / 1000.0
ev["onset"] = ev["onset"] - T0
if DUR is not None:
    ev["duration"] = float(DUR)
ev = ev[(ev["onset"] >= 0) & (ev["onset"] < LEN)]

if MODE == "task":
    dropped = ev["trial_type"].isin(REST)
    print("events modelled as task: %d | left in the baseline (%s): %d"
          % (int((~dropped).sum()), ",".join(REST), int(dropped.sum())))
    ev = ev[~dropped].assign(trial_type="task")
    contrast = "task"
else:
    contrast = MODE
    print("trial types kept:", dict(ev["trial_type"].value_counts()))
ev = ev[["onset", "duration", "trial_type"]]
busy = float(np.clip(ev["duration"], 0.5, None).sum()) / LEN
print("the modelled events cover roughly %.0f%% of the scan" % (100 * busy))

# ---- confounds and fit ------------------------------------------------------
cf = pd.read_csv(CONF, sep="\t")
mot = [c for c in ("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z") if c in cf.columns]
use = mot + [c for c in cf.columns if c.startswith("a_comp_cor")]
conf = cf[use].fillna(0.0)

model = FirstLevelModel(t_r=TR, hrf_model="spm", drift_model="cosine", high_pass=1 / 128.0,
                        noise_model="ar1", smoothing_fwhm=FWHM, mask_img=MASK, minimize_memory=False)
model.fit(BOLD, events=ev, confounds=conf)
z = model.compute_contrast(contrast, output_type="z_score")
nib.save(z, REP + "z_" + contrast.replace(" ", "") + ".nii.gz")

# ---- is the activation just task-correlated head motion? --------------------
dm = model.design_matrices_[0]
reg = dm[contrast] if contrast in dm.columns else dm.iloc[:, 0]
print("\ntask regressor vs motion (|r| above ~0.3 means the task and the head move together):")
worst = 0.0
for c in mot + (["framewise_displacement"] if "framewise_displacement" in cf.columns else []):
    v = pd.to_numeric(cf[c], errors="coerce").fillna(0.0).values[:len(reg)]
    r = float(np.corrcoef(reg.values, v)[0, 1]) if np.std(v) > 0 else 0.0
    worst = max(worst, abs(r))
    print("   %-24s r = %+.3f" % (c, r))
print("   largest |r| = %.3f -> %s" % (worst, "fine" if worst < 0.3 else "CAUTION: check this"))

# ---- numbers and pictures ---------------------------------------------------
zd = np.nan_to_num(z.get_fdata())
inmask = nib.load(MASK).get_fdata() > 0
hit = (zd > THR) & inmask
print("\nmax z %.2f | voxels z>%.1f: %d of %d" % (zd[inmask].max(), THR, int(hit.sum()), int(inmask.sum())))
if os.path.exists(LAB):
    lab = np.rint(nib.load(LAB).get_fdata()).astype(int)
    for name, ids in (("cerebellar GM", (8, 47)), ("cerebellar WM", (7, 46))):
        sel = np.isin(lab, ids) & inmask
        k = int((hit & sel).sum())
        print("   %-14s %6d active of %6d voxels (%.1f%%)" % (name, k, int(sel.sum()),
              100.0 * k / max(int(sel.sum()), 1)))

tab = get_clusters_table(z, stat_threshold=THR, cluster_threshold=10)
tab.to_csv(REP + "clusters.csv", index=False)
print("\ntop clusters:\n", tab.head(10).to_string())

plot_design_matrix(dm)
plt.savefig(REP + "design_matrix.png", dpi=120, bbox_inches="tight")
plt.close("all")

bg = BG if os.path.exists(BG) else None
for axis in ("z", "y", "x"):
    plot_stat_map(z, bg_img=bg, threshold=THR, display_mode=axis, cut_coords=7,
                  black_bg=False, title="%s  %s  z>%.1f" % (TASK, contrast, THR))
    plt.savefig(REP + "map_" + axis + ".png", dpi=120, bbox_inches="tight")
    plt.close("all")

# ---- the HTML page ----------------------------------------------------------
html = REP + "glm_report.html"
try:
    from nilearn.reporting import make_glm_report
    make_glm_report(model, contrasts=contrast, bg_img=bg, threshold=THR,
                    cluster_threshold=10, title="%s  %s" % (TASK, contrast)).save_as_html(html)
except Exception as exc:
    print("make_glm_report failed (%s); trying model.generate_report" % exc)
    try:
        model.generate_report(contrasts=contrast, bg_img=bg, threshold=THR,
                              cluster_threshold=10, title="%s  %s" % (TASK, contrast)).save_as_html(html)
    except Exception as exc2:
        html = None
        print("no HTML report:", exc2)

print("\nwritten to", REP)
if html:
    print("open", html, "in a browser (or the VS Code preview)")
