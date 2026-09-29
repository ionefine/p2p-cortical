"""
Python port of ed.m -- the Elche Drawings pipeline orchestration layer
(patient-drawing I/O, phosphene simulation/matching against patient
drawings, random-phosphene generation, and real-vs-random statistical
comparison).

Mirrors the CURRENT (bug-fixed) state of ElcheDrawings/ed.m in this repo:
  - the premature plot_corr_histograms call (Main script) and the
    idx(i+4) bug (visualize_combinations) found during review have
    already been fixed in the MATLAB source this was ported from, and
    this Python port reflects that fixed behavior.
  - known-intentional design choices confirmed with the author are kept
    exactly as-is, not "fixed":
      * random phosphenes are NOT rotated during registration
        (simulate_save_rand_drawings forces r=0), unlike real singletons.
      * only the FIRST electrode of a multi-electrode "random combo" is
        actually randomized (combine_random_models / show_best_combos);
        the rest stay real. This is intentional: the question being asked
        is "does adding electrode N+1 beat chance," not "is the whole
        combo no better than chance."
      * combine_sim_draws / combine_random_models always regress against
        patient_img[0] (the first subimage), even when pooling candidate
        columns from every subimage.
  - norm255's known quirk (dividing by vmax instead of vmax-vmin in the
    rectify branch) is preserved as-is; it's a no-op whenever vmin==0,
    which is true for every caller here (imwarp's zero-fill background
    guarantees it).

Python vs. MATLAB struct semantics: MATLAB structs are copy-on-write value
types; Python objects are references. Every place ed.m does something like
`c.e = c_orig.e(eIdx)` to get a scratch single-electrode copy, this port
uses `copy.deepcopy` on the electrode struct so later in-place mutation
(setting .ef, .rfmap, etc.) can't leak back into the shared original. Large
read-only arrays on c/v (X, Y, ORmap, cropPix, ...) are intentionally
shared by reference rather than deep-copied -- they're never mutated after
generate_corticalmap, so sharing them is behaviorally equivalent to
MATLAB's copy-on-write and avoids needless duplication of large arrays.

File I/O: MATLAB .mat files are replaced with Python pickle (.pkl) files,
using the same directory layout (models/, combos/, random_models/,
figures/) and filename convention (vbl.fileidstr, debug tag) as the MATLAB
version. The two pipelines' output files are therefore NOT directly
interchangeable, but the naming *logic* is preserved 1:1.

SSIM: implemented as a standard 11x11-Gaussian-windowed structure-only SSIM
(Wang et al. 2004 formula, exponents [0,0,1] i.e. structure term only,
matching MATLAB's `ssim(...,'Exponents',[0 0 1],'DynamicRange',255)` call
sites). This is NOT reverse-engineered from MATLAB's internal ssim()
implementation bit-for-bit; that's fine because, in both the MATLAB and
this Python pipeline, the ssim value is only ever logged/saved, never used
to select or rank candidates (correlation via fastcorr drives every
selection decision).
"""

from __future__ import annotations

import concurrent.futures as cf
import copy
import itertools
import os
import pickle
from pathlib import Path
from typing import List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.signal import convolve2d

import p2p_c
from p2p_c import Struct
from registration import (
    find_scale_rotation_ngc,
    resolve_similarity_rotation_ambiguity_ngc,
    imwarp_similarity_fixed,
)
from crop_img import crop_img


# ===========================================================================
# Small helpers
# ===========================================================================
def norm01(img: np.ndarray) -> np.ndarray:
    img = np.asarray(img, dtype=np.float64)
    rng = img.max() - img.min()
    if rng < np.finfo(float).eps:
        return np.zeros_like(img)
    return (img - img.min()) / rng


def norm255(v: np.ndarray, rectify: bool):
    """Port of ed.m's norm255, quirk included (see module docstring)."""
    v = np.asarray(v, dtype=np.float64).ravel()
    vmax = v.max()
    vmin = v.min()
    if rectify:
        v = np.maximum(v - vmin, 0)
        if vmax > 0:
            v = v / vmax * 255
        return np.clip(np.round(v), 0, 255).astype(np.uint8)
    else:
        denom = max(vmax - vmin, np.finfo(float).eps)
        out = np.round(128 * (v - vmin) / denom)
        return np.clip(out, -128, 127).astype(np.int8)


def fastcorr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    a = a - a.mean()
    b = b - b.mean()
    na = np.sqrt(np.sum(a * a))
    nb = np.sqrt(np.sum(b * b))
    if na == 0 or nb == 0:
        return 0.0
    return float(np.sum(a * b) / (na * nb))


def _ssim_structure(img1: np.ndarray, img2: np.ndarray, dynamic_range: float = 255.0) -> float:
    """Structure-only SSIM (see module docstring for why exact MATLAB
    numeric parity isn't required here)."""
    img1 = np.asarray(img1, dtype=np.float64)
    img2 = np.asarray(img2, dtype=np.float64)

    K2 = 0.03
    C3 = (K2 * dynamic_range) ** 2 / 2

    sigma = 1.5
    radius = int(3.5 * sigma + 0.5)
    ax = np.arange(-radius, radius + 1)
    win1d = np.exp(-(ax**2) / (2 * sigma**2))
    win1d /= win1d.sum()
    win = np.outer(win1d, win1d)

    def filt(x):
        return convolve2d(x, win, mode="valid")

    mu1 = filt(img1)
    mu2 = filt(img2)
    sigma1 = np.maximum(filt(img1 * img1) - mu1 * mu1, 0)
    sigma2 = np.maximum(filt(img2 * img2) - mu2 * mu2, 0)
    sigma12 = filt(img1 * img2) - mu1 * mu2

    s_map = (sigma12 + C3) / (np.sqrt(sigma1) * np.sqrt(sigma2) + C3)
    return float(np.mean(s_map))


def ensure_dirs(draw_dir: Path):
    for sub in ("models", "combos", "random_models", "drawings"):
        (draw_dir / sub).mkdir(parents=True, exist_ok=True)


def clean_dir_safe(p: Path):
    if not p.exists():
        return
    for entry in p.iterdir():
        if entry.is_dir():
            continue
        try:
            entry.unlink()
        except OSError:
            pass


def _tag(vbl) -> str:
    return "_debug" if vbl.debugflag else ""


def _save(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def _load(path: Path):
    with open(path, "rb") as f:
        return pickle.load(f)


def safe_rmfield(struct_obj: Struct, names) -> Struct:
    """Port of ed.m's safe_rmfield: drop `names` from `struct_obj` if present,
    no-op for names that aren't there. Used by the Main script to strip heavy
    fields before shipping c/v to parallel random-rep workers."""
    if struct_obj is None:
        return struct_obj
    for name in names:
        if hasattr(struct_obj, name):
            delattr(struct_obj, name)
    return struct_obj


def _struct_copy_without(struct_obj: Struct, exclude=()) -> Struct:
    """Shallow copy of top-level struct fields, excluding `exclude`. Safe
    because the only fields ever mutated downstream through such a copy are
    the excluded ones (typically 'e'); large arrays (X, Y, ORmap, ...) are
    intentionally shared by reference (never mutated in place)."""
    return Struct(**{k: v for k, v in vars(struct_obj).items() if k not in exclude})


# ===========================================================================
# Setup
# ===========================================================================
def setup(loc: str = "auto") -> Struct:
    vbl = Struct()

    if loc not in ("Z", "auto") and Path(loc).is_dir():
        vbl.datadir = Path(loc)
    else:
        vbl.datadir = Path.cwd()

    vbl.n_SimDraw = 24
    vbl.n_DrawCombined = 3
    vbl.debugflag = 0
    vbl.n_Reps = 4 if vbl.debugflag == 1 else 100
    vbl.rectify = 1
    vbl.fileidstr = "_10_2_2025"

    vbl.cleanstartflag = 1
    resp = input("Delete all previous simulations? 1=yes, 0=no ... ")
    vbl.cleanstartflag = int(resp)

    vbl.dirList = sorted(
        [p for p in Path(vbl.datadir).iterdir() if p.is_dir() and "img" in p.name],
        key=lambda p: p.name,
    )
    vbl.n_Drawings = len(vbl.dirList)
    if vbl.n_Drawings == 0:
        raise FileNotFoundError(f'No drawing folders found in {vbl.datadir} (pattern "*img*").')
    if vbl.debugflag:
        vbl.n_Drawings = 10

    return vbl


# ===========================================================================
# Patient drawings I/O
# ===========================================================================
def load_patient_drawings(vbl) -> List[Struct]:
    p_draw = []
    for d in range(vbl.n_Drawings):
        draw_dir = Path(vbl.datadir) / vbl.dirList[d].name
        print(f"loading {draw_dir}")

        if vbl.cleanstartflag:
            ensure_dirs(draw_dir)
            clean_dir_safe(draw_dir / "combos")
            clean_dir_safe(draw_dir / "models")
            clean_dir_safe(draw_dir / "random_models")
            clean_dir_safe(draw_dir)  # remove all files at drawing root (subfolders untouched)
        else:
            ensure_dirs(draw_dir)

        sub_imgs = sorted((draw_dir / "drawings").glob("*edit_crop_*.png"))
        pd = Struct()
        pd.n_SubImages = len(sub_imgs)
        if pd.n_SubImages == 0:
            raise FileNotFoundError(f"No subimages found in {draw_dir}/drawings")

        pd.nToSave = int(np.ceil(vbl.n_SimDraw / pd.n_SubImages))

        pd.patient_img = {}
        pd.size = {}
        pd.sim_img = {}
        pd.corr = {}
        pd.ssim = {}
        pd.radius = {}
        pd.x = {}
        pd.y = {}
        pd.subID = {}
        pd.subimg = {}
        pd.canonicalSize = None
        pd.ref2d = None

        for sd in range(pd.n_SubImages):
            fname = draw_dir / "drawings" / f"{vbl.dirList[d].name}_edit_crop_{sd}.png"
            img = np.array(Image.open(fname))
            if img.ndim == 3:
                img = img.mean(axis=2)
            img = np.clip(np.round(img), 0, 255).astype(np.uint8)

            pd.patient_img[sd] = img
            pd.size[sd] = img.shape

            if sd == 0:
                pd.ref2d = img.shape
                pd.canonicalSize = img.shape
            else:
                if img.shape != pd.canonicalSize:
                    raise ValueError(
                        "load_patient_drawings:subimageSizeMismatch\n"
                        f"Subimage sizes differ within drawing {vbl.dirList[d].name} (d={d}).\n"
                        f"sd=0 size: {pd.canonicalSize}\nsd={sd} size: {img.shape}"
                    )

            K = pd.nToSave
            pd.sim_img[sd] = np.zeros(img.shape + (vbl.n_SimDraw,), dtype=np.uint8)
            pd.corr[sd] = np.full(K, -np.inf)
            pd.ssim[sd] = np.full(K, -np.inf)
            pd.radius[sd] = np.full(K, np.nan)
            pd.x[sd] = np.full(K, np.nan)
            pd.y[sd] = np.full(K, np.nan)
            pd.subID[sd] = np.full(K, np.nan)
            pd.subimg[sd] = np.full((img.size, K), np.nan)

        p_draw.append(pd)
    return p_draw


def crop_drawings(vbl):
    for d in range(vbl.n_Drawings):
        draw_name = vbl.dirList[d].name
        draw_dir = Path(vbl.datadir) / draw_name / "drawings"

        orig_file = draw_dir / f"{draw_name}.png"
        edit_file = draw_dir / f"{draw_name}_edit.png"

        if not orig_file.is_file() or not edit_file.is_file():
            print(f"Missing original or edited image in {draw_dir}")
            continue

        orig_img = np.array(Image.open(orig_file), dtype=np.float64)
        edit_img = np.array(Image.open(edit_file), dtype=np.float64)

        img = orig_img.mean(axis=2) if orig_img.ndim == 3 else orig_img
        sz = img.shape

        row_mid = round(sz[0] / 2)
        col_mid = round(sz[1] / 2)

        cidx = np.flatnonzero(img[row_mid, :] == 0)
        ridx = np.flatnonzero(img[:, col_mid] == 0)

        if len(cidx) < 2 or len(ridx) < 2:
            print(f"Failed to auto-crop {draw_name}; saving originals as crops.")
            orig_img_tmp = orig_img
            edit_img_tmp = edit_img
        else:
            r1, c1, r2, c2 = ridx[0] + 1, cidx[0] + 1, ridx[-1] - 1, cidx[-1] - 1
            orig_img_tmp = orig_img[r1 : r2 + 1, c1 : c2 + 1, ...]
            edit_img_tmp = edit_img[r1 : r2 + 1, c1 : c2 + 1, ...]

        Image.fromarray((norm01(orig_img_tmp) * 255).astype(np.uint8)).save(draw_dir / f"{draw_name}_crop.png")
        Image.fromarray((norm01(edit_img_tmp) * 255).astype(np.uint8)).save(
            draw_dir / f"{draw_name}_edit_crop_0.png"
        )


# ===========================================================================
# Cortical/visual model
# ===========================================================================
def define_cortical_model(vbl):
    print("Defining cortical/visual model...")

    tp = p2p_c.define_temporalparameters()
    trl = Struct(freq=float("nan"))
    trl = p2p_c.define_trial(tp, trl)

    c = Struct()
    v = Struct()
    if vbl.debugflag:
        v.pixperdeg = 4
        c.pixpermm = 4
        n_electrodes = 20
    else:
        v.pixperdeg = 10
        c.pixpermm = 10
        n_electrodes = 1000

    c.cortexHeight = [-30, 30]
    c.cortexLength = [-55, 55]
    c.onoff_ratio = 0.9

    v.visfieldHeight = [-30, 30]
    v.visfieldWidth = [-30, 30]

    v = p2p_c.define_visualmap(v)
    c = p2p_c.define_cortex(c)
    c, v = p2p_c.generate_corticalmap(c, v)

    theta = 2 * np.pi * np.random.rand(n_electrodes)
    r = np.exp(4 * np.random.rand(n_electrodes))
    r = r * 3 / np.max(r)
    ecc = np.exp(np.log(15) * np.random.rand(n_electrodes))
    c.I_k = 1000

    x = ecc * np.cos(theta)
    y = ecc * np.sin(theta)
    iperm = np.random.permutation(n_electrodes)

    v.e = [Struct(x=float(x[iperm[i]]), y=float(y[iperm[i]])) for i in range(n_electrodes)]
    c.e = [Struct(radius=float(r[iperm[i]])) for i in range(n_electrodes)]

    c, v = p2p_c.define_electrodes(c, v)
    print(f"Defined {n_electrodes} electrodes.")

    return c, v, trl, tp


# ===========================================================================
# Simulation of singletons
# ===========================================================================
_worker_state: dict = {}


def _init_electrode_worker(c_orig, v_orig, trl, tp, targets, vbl):
    """ProcessPoolExecutor initializer: stashes the (large, read-only)
    shared inputs in this worker process's globals ONCE, instead of
    re-pickling/re-sending them on every task submission."""
    _worker_state["c_orig"] = c_orig
    _worker_state["v_orig"] = v_orig
    _worker_state["trl"] = trl
    _worker_state["tp"] = tp
    _worker_state["targets"] = targets
    _worker_state["vbl"] = vbl


def _compute_electrode_candidates(e_idx: int):
    """Phase-1 worker: everything about a single electrode that does NOT
    depend on shared/mutable pool state -- phosphene generation (the
    profiled bottleneck) plus registration against every (drawing,
    subimage) target. Fully independent across electrodes, so safe to run
    in parallel/in any order; the electrode-order-dependent pool-selection
    logic lives in simulate_drawings's phase-2 merge loop instead."""
    c_orig = _worker_state["c_orig"]
    v_orig = _worker_state["v_orig"]
    trl = _worker_state["trl"]
    tp = _worker_state["tp"]
    targets = _worker_state["targets"]
    vbl = _worker_state["vbl"]

    c = _struct_copy_without(c_orig, exclude=("e",))
    v = _struct_copy_without(v_orig, exclude=("e",))
    c.e = [copy.deepcopy(c_orig.e[e_idx])]
    v.e = [copy.deepcopy(v_orig.e[e_idx])]

    c = p2p_c.generate_ef(c)
    v, c = p2p_c.generate_corticalelectricalresponse(c, v)
    img = generate_phosphene(v, tp, trl, vbl)

    candidates = {}
    for d, sd, target, target_size, ref2d in targets:
        s, r = find_scale_rotation_ngc(img.astype(np.float32), target.astype(np.float32))
        A, peakcorr = resolve_similarity_rotation_ambiguity_ngc(
            img.astype(np.float32), target.astype(np.float32), s, r
        )
        img_aligned = imwarp_similarity_fixed(img, A, ref2d, fill_value=0)

        if img_aligned.shape != target_size:
            raise ValueError(
                "simulate_drawings:sizeMismatch "
                f"Warp output size != target size (e_idx={e_idx}, d={d}, sd={sd})."
            )

        pix_count = target.size
        scaled_vec = norm255(img_aligned.ravel(), vbl.rectify).astype(np.float64)
        if scaled_vec.size != pix_count:
            raise ValueError("simulate_drawings:vectorLengthMismatch Vector length mismatch after warp.")

        new_corr = fastcorr(target.ravel(), img_aligned.ravel())
        new_ssim = _ssim_structure(img_aligned.reshape(target_size), target, 255)

        candidates[(d, sd)] = {
            "peakcorr": peakcorr,
            "scaled_vec": scaled_vec,
            "new_corr": new_corr,
            "new_ssim": new_ssim,
            "target_size": target_size,
        }

    return e_idx, candidates, c.e[0].radius, v.e[0].x, v.e[0].y


def simulate_drawings(c_orig, v_orig, trl, tp, p_draw, vbl, max_workers=None):
    """Parallelized across electrodes (each electrode's phosphene
    generation -- p2p_c.generate_corticalelectricalresponse, profiled at
    ~130s/electrode at production resolution, i.e. tens of hours for a
    1000-electrode run done serially -- and registration against every
    drawing/subimage is fully independent of every other electrode).

    This is split into two phases to guarantee byte-identical output to a
    fully-serial version:
      1. `_compute_electrode_candidates` runs in a ProcessPoolExecutor,
         one task per electrode, computing everything that doesn't touch
         shared mutable state.
      2. A plain, strictly-sequential loop over electrodes in original
         (0..n_elect-1) order applies the pool-selection logic (the "quick
         gate" and diversity check, which compare against whatever's
         already been saved, so are genuinely order-dependent) to the
         phase-1 results.

    Must be called from code guarded by `if __name__ == "__main__":` when
    using the default process-pool backend (already true of
    run_elche_drawings_main.py's main()).
    """
    corr_thresh = 0.90
    n_elect = len(c_orig.e)
    if n_elect == 0:
        return p_draw

    targets = []
    for d in range(vbl.n_Drawings):
        pd = p_draw[d]
        for sd in range(pd.n_SubImages):
            targets.append((d, sd, pd.patient_img[sd], pd.patient_img[sd].shape, pd.ref2d))

    if max_workers is None:
        max_workers = os.cpu_count() or 1
    max_workers = max(1, min(max_workers, n_elect))

    results_by_electrode = {}
    with cf.ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_init_electrode_worker,
        initargs=(c_orig, v_orig, trl, tp, targets, vbl),
    ) as executor:
        futures = [executor.submit(_compute_electrode_candidates, e_idx) for e_idx in range(n_elect)]
        n_done = 0
        for f in cf.as_completed(futures):
            e_idx, candidates, radius, elec_x, elec_y = f.result()
            results_by_electrode[e_idx] = (candidates, radius, elec_x, elec_y)
            n_done += 1
            print(f"Electrode {e_idx + 1} / {n_elect} computed ({n_done}/{n_elect} done)")

    # Phase 2: sequential merge in original electrode order -- preserves
    # exact selection semantics of the old fully-serial version.
    for e_idx in range(n_elect):
        candidates, radius, elec_x, elec_y = results_by_electrode[e_idx]

        for d in range(vbl.n_Drawings):
            pd = p_draw[d]
            for sd in range(pd.n_SubImages):
                cand = candidates[(d, sd)]
                target_size = cand["target_size"]
                peakcorr = cand["peakcorr"]

                if peakcorr <= np.min(pd.corr[sd]):
                    continue

                pix_count = int(np.prod(target_size))
                if pd.subimg[sd].shape[0] != pix_count:
                    raise ValueError("simulate_drawings:poolSizeMismatch Pool pixel dimension mismatch BEFORE write.")

                scaled_vec = cand["scaled_vec"]
                new_corr = cand["new_corr"]
                new_ssim = cand["new_ssim"]

                pool = pd.subimg[sd]
                K = pool.shape[1]

                if np.all(~np.isfinite(pd.corr[sd])) or np.all(np.isnan(pool)):
                    max_pool_corr = 0.0
                    most_similar_idx = 0
                else:
                    pc = np.full(K, -np.inf)
                    for k in range(K):
                        if np.all(np.isfinite(pool[:, k])):
                            pc[k] = fastcorr(pool[:, k], scaled_vec)
                    most_similar_idx = int(np.argmax(pc))
                    max_pool_corr = pc[most_similar_idx]
                    if not np.isfinite(max_pool_corr):
                        max_pool_corr = 0.0
                        most_similar_idx = 0

                if max_pool_corr > corr_thresh:
                    ridx = most_similar_idx
                else:
                    ridx = int(np.argmin(pd.corr[sd]))

                if new_corr > pd.corr[sd][ridx]:
                    pool[:, ridx] = scaled_vec
                    pd.radius[sd][ridx] = radius
                    pd.x[sd][ridx] = elec_x
                    pd.y[sd][ridx] = elec_y
                    pd.corr[sd][ridx] = new_corr
                    pd.ssim[sd][ridx] = new_ssim
                    pd.subID[sd][ridx] = sd
                    pd.size[sd] = target_size

    return p_draw


def generate_phosphene(v: Struct, tp: Struct, trl: Struct, vbl: Struct) -> np.ndarray:
    trl_array = p2p_c.generate_phosphene(v, tp, trl)
    img = np.max(trl_array[0].max_phosphene, axis=2)
    img, *_ = crop_img(img, 20)

    if img.size == 0:
        # MATLAB: max(abs(img(:))) on an empty img returns [], `if [] < eps`
        # is falsy (skips the m=1 guard), and img./[] stays empty -- so a
        # fully degenerate (uniform/no-response) electrode propagates as an
        # empty array rather than raising. Mirror that here instead of
        # crashing on np.max([]).
        return img.astype(np.uint8)

    m = np.max(np.abs(img))
    if m < np.finfo(float).eps:
        m = 1.0
    img = img / m

    if vbl.rectify:
        img = np.where(img < 0, 0, img)
        img = np.clip(np.round(255 * img), 0, 255).astype(np.uint8)
    else:
        img = np.clip(np.round(127 * (img + 0.5)), 0, 255).astype(np.uint8)

    return img


def save_simulated_drawings(p_draw, vbl):
    for d in range(vbl.n_Drawings):
        pd = p_draw[d]
        sim_draw = Struct()
        sim_draw.n_SubImages = pd.n_SubImages
        sim_draw.nToSave = pd.nToSave
        sim_draw.patient_img = pd.patient_img
        sim_draw.ref2d = pd.ref2d

        sim_draw.subimg = {}
        sim_draw.radius = {}
        sim_draw.x = {}
        sim_draw.y = {}
        sim_draw.corr = {}
        sim_draw.ssim = {}
        sim_draw.subID = {}
        sim_draw.size = {}

        for sd in range(sim_draw.n_SubImages):
            order = np.argsort(-pd.corr[sd])
            K = min(sim_draw.nToSave, len(order))

            sim_draw.subimg[sd] = np.zeros((pd.subimg[sd].shape[0], K))
            sim_draw.radius[sd] = np.zeros(K)
            sim_draw.x[sd] = np.zeros(K)
            sim_draw.y[sd] = np.zeros(K)
            sim_draw.corr[sd] = np.zeros(K)
            sim_draw.ssim[sd] = np.zeros(K)

            for i in range(K):
                ii = order[i]
                sim_draw.subimg[sd][:, i] = pd.subimg[sd][:, ii]
                sim_draw.radius[sd][i] = pd.radius[sd][ii]
                sim_draw.x[sd][i] = pd.x[sd][ii]
                sim_draw.y[sd][i] = pd.y[sd][ii]
                sim_draw.corr[sd][i] = pd.corr[sd][ii]
                sim_draw.ssim[sd][i] = pd.ssim[sd][ii]

            sim_draw.subID[sd] = sd
            sim_draw.size[sd] = pd.size[sd]

        out_dir = Path(vbl.datadir) / vbl.dirList[d].name / "models"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_name = f"{vbl.dirList[d].name}{vbl.fileidstr}{_tag(vbl)}.pkl"
        _save(sim_draw, out_dir / out_name)


# ===========================================================================
# Random phosphenes
# ===========================================================================
def simulate_save_rand_drawings(rep, c, v, trl, tp, vbl):
    print(f"Random rep {rep + 1} / {vbl.n_Reps}")

    c2, v2 = p2p_c.generate_corticalmap(c, v)
    tag = _tag(vbl)

    for d in range(vbl.n_Drawings):
        draw_name = vbl.dirList[d].name
        draw_dir = Path(vbl.datadir) / draw_name
        models_dir = draw_dir / "models"
        rand_dir = draw_dir / "random_models"
        rand_dir.mkdir(parents=True, exist_ok=True)

        out_file = rand_dir / f"{draw_name}{vbl.fileidstr}_{rep}{tag}.pkl"
        if out_file.exists():
            continue  # already saved for this drawing+rep

        sim_file = models_dir / f"{draw_name}{vbl.fileidstr}{tag}.pkl"
        if not sim_file.is_file():
            raise FileNotFoundError(f"simulate_save_rand_drawings:missingSimDraw Missing sim_draw file: {sim_file}")

        sim_draw = _load(sim_file)

        rand_draw = Struct()
        rand_draw.subimg = {}
        rand_draw.corr = {}
        rand_draw.ssim = {}

        for sd in range(sim_draw.n_SubImages):
            print(f"{draw_name} subdrawing {sd}")
            target = sim_draw.patient_img[sd]
            target_size = target.shape
            pix_count = target.size
            K = sim_draw.nToSave

            rand_draw.subimg[sd] = np.full((pix_count, K), np.nan)
            rand_draw.corr[sd] = np.full(K, np.nan)
            rand_draw.ssim[sd] = np.full(K, np.nan)

            for i in range(K):
                c2e = _struct_copy_without(c2, exclude=("e",))
                v2e = _struct_copy_without(v2, exclude=("e",))
                c2e.e = [Struct(radius=float(sim_draw.radius[sd][i]))]
                v2e.e = [Struct(x=float(sim_draw.x[sd][i]), y=float(sim_draw.y[sd][i]))]

                c2e, _ = p2p_c.define_electrodes(c2e, v2e)
                c2e = p2p_c.generate_ef(c2e)
                v2e, c2e = p2p_c.generate_corticalelectricalresponse(c2e, v2e)

                img = generate_phosphene(v2e, tp, trl, vbl)

                # Align to THIS subimage sd
                s, r = find_scale_rotation_ngc(img.astype(np.float32), target.astype(np.float32))

                # Random phosphenes are intentionally NOT rotated (confirmed
                # intentional, unlike real singletons in simulate_drawings).
                r = 0.0

                A, _peak = resolve_similarity_rotation_ambiguity_ngc(
                    img.astype(np.float32), target.astype(np.float32), s, r
                )
                img_aligned = imwarp_similarity_fixed(img, A, target_size, fill_value=0)

                if img_aligned.shape != target_size:
                    raise ValueError("simulate_save_rand_drawings:sizeMismatch Warped random size mismatch.")

                scaled_vec = norm255(img_aligned.ravel(), vbl.rectify).astype(np.float64)
                if scaled_vec.size != pix_count:
                    raise ValueError(
                        "simulate_save_rand_drawings:vectorLengthMismatch Vector length mismatch after warp."
                    )

                rand_draw.subimg[sd][:, i] = scaled_vec
                rand_draw.corr[sd][i] = fastcorr(target.ravel(), img_aligned.ravel())
                rand_draw.ssim[sd][i] = _ssim_structure(img_aligned.reshape(target_size), target, 255)

        _save(rand_draw, out_file)


# ===========================================================================
# Combine real and random
# ===========================================================================
def _enumerate_combo_masks(n_cols: int, max_size: int):
    """All boolean masks of size n_cols with 1..max_size bits set -- the
    same set MATLAB's dec2bin(1:(2^nCols-1)) + sum<=max_size filter
    produces, generated directly instead of via a 2^nCols brute-force scan
    (an efficiency change only; the resulting combo set is identical)."""
    masks = []
    for k in range(1, max_size + 1):
        for idxs in itertools.combinations(range(n_cols), k):
            mask = np.zeros(n_cols, dtype=bool)
            mask[list(idxs)] = True
            masks.append(mask)
    return masks


def combine_sim_draws(vbl):
    for d in range(vbl.n_Drawings):
        draw_dir = Path(vbl.datadir) / vbl.dirList[d].name
        models_dir = draw_dir / "models"
        combos_dir = draw_dir / "combos"
        combos_dir.mkdir(parents=True, exist_ok=True)

        print(f"Combining real singletons: {d + 1} / {vbl.n_Drawings}")

        tag = _tag(vbl)
        in_file = models_dir / f"{vbl.dirList[d].name}{vbl.fileidstr}{tag}.pkl"
        sim_draw = _load(in_file)

        columns = [
            sim_draw.subimg[sd][:, ii].astype(np.float64)
            for sd in range(sim_draw.n_SubImages)
            for ii in range(sim_draw.nToSave)
        ]
        tmp_subimg = np.column_stack(columns)
        n_cols = tmp_subimg.shape[1]

        combo_masks = _enumerate_combo_masks(n_cols, vbl.n_DrawCombined)
        nC = len(combo_masks)

        nimg = np.zeros(nC)
        corr_v = np.zeros(nC)
        ssim_v = np.zeros(nC)
        cmbx = np.zeros((nC, n_cols), dtype=bool)

        target_img = sim_draw.patient_img[0].astype(np.float64)
        target = target_img.ravel()

        for ci, mask in enumerate(combo_masks):
            X = tmp_subimg[:, mask] / 255.0
            B, *_ = np.linalg.lstsq(X, target, rcond=None)
            recon = X @ B
            nimg[ci] = mask.sum()
            corr_v[ci] = fastcorr(recon, target)
            ssim_v[ci] = _ssim_structure(recon.reshape(target_img.shape), target_img, 255)
            cmbx[ci, :] = mask

        order = np.argsort(-corr_v)
        combo = Struct()
        combo.nimg = nimg[order]
        combo.corr_val = corr_v[order]
        combo.ssim_val = ssim_v[order]
        combo.cmbx = cmbx[order, :]
        # Note: combo.subimg omitted to reduce disk footprint (matches MATLAB ed.m)

        out_file = combos_dir / f"{vbl.dirList[d].name}{vbl.fileidstr}{tag}.pkl"
        _save(combo, out_file)


def combine_random_models(vbl):
    for d in range(vbl.n_Drawings):
        draw_name = vbl.dirList[d].name
        draw_dir = Path(vbl.datadir) / draw_name
        combos_dir = draw_dir / "combos"
        models_dir = draw_dir / "models"

        print(f"Combining random models: {d + 1} / {vbl.n_Drawings}")

        tag = _tag(vbl)
        combo_file = combos_dir / f"{draw_name}{vbl.fileidstr}{tag}.pkl"
        sim_file = models_dir / f"{draw_name}{vbl.fileidstr}{tag}.pkl"

        combo = _load(combo_file)
        sim_draw = _load(sim_file)

        columns = [
            sim_draw.subimg[sd][:, ii].astype(np.float64)
            for sd in range(sim_draw.n_SubImages)
            for ii in range(sim_draw.nToSave)
        ]
        tmp_subimg = np.column_stack(columns)

        nComb = combo.cmbx.shape[0]
        R = vbl.n_Reps
        rand_combo = Struct()
        rand_combo.nimg = np.full((R, nComb), np.nan)
        rand_combo.corr_val = np.full((R, nComb), np.nan)
        rand_combo.ssim_val = np.full((R, nComb), np.nan)
        rand_combo.id = [None] * R
        rand_combo.perms = [None] * R
        rand_combo.cmbx = combo.cmbx

        target_img = sim_draw.patient_img[0].astype(np.float64)
        target = target_img.ravel()

        for rep in range(R):
            tag_r = f"_{rep}_debug" if vbl.debugflag else f"_{rep}"
            rand_file = draw_dir / "random_models" / f"{draw_name}{vbl.fileidstr}{tag_r}.pkl"

            if not rand_file.is_file():
                print(f"Missing random file for rep {rep} in {draw_name}")
                continue
            rand_draw = _load(rand_file)

            randimg_cols = []
            perms_rep = {}
            for sd in range(len(rand_draw.subimg)):
                nK = len(rand_draw.corr[sd])
                shuf_ind = np.random.permutation(nK)
                perms_rep[sd] = shuf_ind
                for ii in shuf_ind:
                    randimg_cols.append(rand_draw.subimg[sd][:, ii].astype(np.float64))
            randimg_local = np.column_stack(randimg_cols)

            local_corr = np.full(nComb, np.nan)
            local_ssim = np.full(nComb, np.nan)
            local_nimg = np.full(nComb, np.nan)

            for ci in range(nComb):
                mask = combo.cmbx[ci, :]
                if not mask.any():
                    continue
                X = (tmp_subimg[:, mask] / 255.0).copy()
                Xr = randimg_local[:, mask] / 255.0

                # Overwrite first column with random (intentional; see module docstring)
                X[:, 0] = Xr[:, 0]

                B, *_ = np.linalg.lstsq(X, target, rcond=None)
                recon = X @ B

                local_nimg[ci] = mask.sum()
                local_corr[ci] = fastcorr(recon, target)
                local_ssim[ci] = _ssim_structure(recon.reshape(target_img.shape), target_img, 255)

            order = np.argsort(-local_corr)
            rand_combo.nimg[rep, :] = local_nimg[order]
            rand_combo.corr_val[rep, :] = local_corr[order]
            rand_combo.ssim_val[rep, :] = local_ssim[order]
            rand_combo.id[rep] = order
            rand_combo.perms[rep] = perms_rep

        out_file = combos_dir / f"{draw_name}{vbl.fileidstr}_rand{tag}.pkl"
        _save(rand_combo, out_file)


# ===========================================================================
# Visualization
# ===========================================================================
def visualize_singletons(vbl, d):
    tag = _tag(vbl)
    in_file = Path(vbl.datadir) / vbl.dirList[d].name / "models" / f"{vbl.dirList[d].name}{vbl.fileidstr}{tag}.pkl"
    if not in_file.is_file():
        raise FileNotFoundError(f"Missing file: {in_file}")
    sim_draw = _load(in_file)

    fig, axes = plt.subplots(2, 2)
    for i, ax in enumerate(axes.flat):
        if i < sim_draw.n_SubImages:
            ax.imshow(sim_draw.patient_img[i], cmap="gray")
        ax.axis("off")
    fig.suptitle(vbl.dirList[d].name)

    for s in range(sim_draw.n_SubImages):
        fig2, axes2 = plt.subplots(3, 2)
        K = min(6, sim_draw.subimg[s].shape[1])
        for i, ax in enumerate(axes2.flat):
            if i < K:
                ax.imshow(sim_draw.subimg[s][:, i].reshape(sim_draw.size[s]), cmap="gray")
                ax.set_title(f"corr={sim_draw.corr[s][i]:.3f}")
            ax.axis("off")
        fig2.suptitle(f"{vbl.dirList[d].name} subimage {s}")


def visualize_combinations(vbl, d):
    tag = _tag(vbl)
    draw_dir = Path(vbl.datadir) / vbl.dirList[d].name
    sim_file = draw_dir / "models" / f"{vbl.dirList[d].name}{vbl.fileidstr}{tag}.pkl"
    cmb_file = draw_dir / "combos" / f"{vbl.dirList[d].name}{vbl.fileidstr}{tag}.pkl"

    if not sim_file.is_file() or not cmb_file.is_file():
        print(f"Missing models or combos for {vbl.dirList[d].name}")
        return

    sim_draw = _load(sim_file)
    combo = _load(cmb_file)

    fig, axes = plt.subplots(2, 2)
    for i, ax in enumerate(axes.flat):
        if i < sim_draw.n_SubImages:
            ax.imshow(sim_draw.patient_img[i], cmap="gray")
        ax.axis("off")
    fig.suptitle(vbl.dirList[d].name)

    columns = [
        sim_draw.subimg[sd][:, ii].astype(np.float64)
        for sd in range(sim_draw.n_SubImages)
        for ii in range(sim_draw.nToSave)
    ]
    tmp_subimg = np.column_stack(columns)
    target_img = sim_draw.patient_img[0].astype(np.float64)
    target = target_img.ravel()

    for ni in range(1, 4):
        fig2, axes2 = plt.subplots(2, 3)
        idx = np.flatnonzero(combo.nimg == ni)
        K = min(len(idx), 6)
        for i, ax in enumerate(axes2.flat):
            if i < K:
                mask = combo.cmbx[idx[i], :]
                X = tmp_subimg[:, mask] / 255.0
                B, *_ = np.linalg.lstsq(X, target, rcond=None)
                recon = (X @ B).reshape(target_img.shape)
                ax.imshow(recon, cmap="gray")
                ax.set_title(f"corr={combo.corr_val[idx[i]]:.3f}")
            ax.axis("off")
        fig2.suptitle(f"Best {ni}")


def show_best_combos(vbl, d):
    draw_name = vbl.dirList[d].name
    tag = _tag(vbl)
    draw_dir = Path(vbl.datadir) / draw_name

    combo = _load(draw_dir / "combos" / f"{draw_name}{vbl.fileidstr}{tag}.pkl")
    rand_combo = _load(draw_dir / "combos" / f"{draw_name}{vbl.fileidstr}_rand{tag}.pkl")
    sim_draw = _load(draw_dir / "models" / f"{draw_name}{vbl.fileidstr}{tag}.pkl")

    target_img = sim_draw.patient_img[0].astype(np.float64)
    target_vec = target_img.ravel()

    columns = [
        sim_draw.subimg[sd][:, ii].astype(np.float64)
        for sd in range(sim_draw.n_SubImages)
        for ii in range(sim_draw.nToSave)
    ]
    tmp_subimg = np.column_stack(columns)

    fig, axes = plt.subplots(3, 3, figsize=(9, 9))
    fig.suptitle(f"Drawing {d}: {draw_name}")

    for row, nimg in enumerate((1, 2, 3)):
        recon_real = np.full(target_img.shape, np.nan)
        real_corr_val = np.nan
        real_corr_recomputed = np.nan

        mask_real = combo.nimg == nimg
        if mask_real.any():
            real_idxs = np.flatnonzero(mask_real)
            real_idx = real_idxs[np.argmax(combo.corr_val[mask_real])]
            mask = combo.cmbx[real_idx, :]
            X = tmp_subimg[:, mask] / 255.0
            B, *_ = np.linalg.lstsq(X, target_vec, rcond=None)
            recon_real = (X @ B).reshape(target_img.shape)
            real_corr_val = combo.corr_val[real_idx]
            real_corr_recomputed = fastcorr(recon_real.ravel(), target_vec)

        best_corr = -np.inf
        best_rep = None
        best_pos = None
        for rep in range(vbl.n_Reps):
            this_nimg = rand_combo.nimg[rep, :]
            this_corr = rand_combo.corr_val[rep, :]
            mask = this_nimg == nimg
            if mask.any():
                idxs = np.flatnonzero(mask)
                relpos = np.argmax(this_corr[mask])
                cmax = this_corr[mask][relpos]
                if cmax > best_corr:
                    best_corr = cmax
                    best_rep = rep
                    best_pos = idxs[relpos]

        recon_rand = np.full(target_img.shape, np.nan)
        rand_corr_val = np.nan
        rand_corr_recomputed = np.nan
        if best_rep is not None:
            rand_corr_val = rand_combo.corr_val[best_rep, best_pos]
            rand_c = rand_combo.id[best_rep][best_pos]
            mask = combo.cmbx[rand_c, :]

            tag_r = f"_{best_rep}_debug" if vbl.debugflag else f"_{best_rep}"
            rand_draw_file = draw_dir / "random_models" / f"{draw_name}{vbl.fileidstr}{tag_r}.pkl"
            rand_draw = _load(rand_draw_file)

            randimg_cols = []
            for sd in range(len(rand_draw.subimg)):
                shuf_ind = rand_combo.perms[best_rep][sd]
                for ii in shuf_ind:
                    randimg_cols.append(rand_draw.subimg[sd][:, ii].astype(np.float64))
            randimg_local = np.column_stack(randimg_cols)

            X = (tmp_subimg[:, mask] / 255.0).copy()
            Xr = randimg_local[:, mask] / 255.0
            if mask.any():
                X[:, 0] = Xr[:, 0]
            Br, *_ = np.linalg.lstsq(X, target_vec, rcond=None)
            recon_rand = (X @ Br).reshape(target_img.shape)
            rand_corr_recomputed = fastcorr(recon_rand.ravel(), target_vec)

        axes[row, 0].imshow(target_img, cmap="gray")
        axes[row, 0].axis("off")
        axes[row, 0].set_ylabel(f"nimg = {nimg}")
        if row == 0:
            axes[row, 0].set_title("Target")

        axes[row, 1].imshow(recon_real, cmap="gray")
        axes[row, 1].axis("off")
        axes[row, 1].set_title(f"Real (saved={real_corr_val:.3f}, rec={real_corr_recomputed:.3f})")

        axes[row, 2].imshow(recon_rand, cmap="gray")
        axes[row, 2].axis("off")
        axes[row, 2].set_title(f"Rand (saved={rand_corr_val:.3f}, rec={rand_corr_recomputed:.3f})")

    outdir = Path(vbl.datadir) / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"Drawing_{draw_name}.pdf")
    plt.close(fig)


def plot_corr_histograms(vbl, d):
    draw_name = vbl.dirList[d].name
    tag = _tag(vbl)
    draw_dir = Path(vbl.datadir) / draw_name

    combo = _load(draw_dir / "combos" / f"{draw_name}{vbl.fileidstr}{tag}.pkl")
    rand_combo = _load(draw_dir / "combos" / f"{draw_name}{vbl.fileidstr}_rand{tag}.pkl")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle(f"Correlation Histograms: {draw_name}")

    for i, nimg in enumerate((1, 2, 3)):
        mask_real = combo.nimg == nimg
        best_real_corr = np.max(combo.corr_val[mask_real]) if mask_real.any() else np.nan

        best_rand_corrs = np.full(vbl.n_Reps, np.nan)
        for rep in range(vbl.n_Reps):
            mask = rand_combo.nimg[rep, :] == nimg
            if mask.any():
                best_rand_corrs[rep] = np.max(rand_combo.corr_val[rep, mask])

        ax = axes[i]
        ax.hist(best_rand_corrs[~np.isnan(best_rand_corrs)], color=(0.3, 0.3, 0.8), edgecolor="k", alpha=0.6)
        if not np.isnan(best_real_corr):
            ax.axvline(best_real_corr, color="r", linewidth=2)
        ax.set_xlabel("Correlation")
        ax.set_ylabel("# Reps")
        ax.set_title(f"nimg = {nimg}\nReal={best_real_corr:.3f}")

    outdir = Path(vbl.datadir) / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"Drawing_{draw_name}_histogram.pdf")
    plt.close(fig)


def compare_real_vs_rand_stats(vbl, d):
    draw_name = vbl.dirList[d].name
    tag = _tag(vbl)
    draw_dir = Path(vbl.datadir) / draw_name

    combo = _load(draw_dir / "combos" / f"{draw_name}{vbl.fileidstr}{tag}.pkl")
    rand_combo = _load(draw_dir / "combos" / f"{draw_name}{vbl.fileidstr}_rand{tag}.pkl")

    for nimg in (1, 2, 3):
        mask_real = combo.nimg == nimg
        best_real_corr = np.max(combo.corr_val[mask_real]) if mask_real.any() else np.nan

        best_rand_per_rep = []
        for rep in range(vbl.n_Reps):
            mask = rand_combo.nimg[rep, :] == nimg
            if mask.any():
                best_rand_per_rep.append(np.max(rand_combo.corr_val[rep, mask]))
        best_rand_per_rep = np.array(best_rand_per_rep)

        # MATLAB's prctile([], ...) degrades gracefully to NaN (nimg=3 has no
        # combos at all when vbl.n_DrawCombined < 3, e.g. in debug/small
        # runs); np.percentile raises on an empty array, so guard explicitly
        # to preserve that graceful MATLAB behavior instead of crashing.
        if best_rand_per_rep.size:
            top5 = np.percentile(best_rand_per_rep, 95)
            top1 = np.percentile(best_rand_per_rep, 99)
        else:
            top5 = np.nan
            top1 = np.nan

        print(f"Image {draw_name} | nimg={nimg}")
        print(f"Best REAL corr: {best_real_corr:.4f}")
        print(f"Random top 5% threshold: {top5:.4f}")
        print(f"Random top 1% threshold: {top1:.4f}")

        if best_real_corr >= top1:
            print("Real is within TOP 1% of randoms.")
        elif best_real_corr >= top5:
            print("Real is within TOP 5% of randoms.")
        else:
            print("Real is BELOW top 5% of randoms.")

        if best_rand_per_rep.size:
            percentile = float(np.mean(best_rand_per_rep <= best_real_corr) * 100)
        else:
            percentile = np.nan
        print(f"Real correlation is at the {percentile:.2f} percentile\n")
