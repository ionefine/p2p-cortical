"""
Python port of Simulate_Elche_Drawings_MainFast_5_3_2026.m -- the entry-point
script that drives the full Elche Drawings pipeline end to end:
load/crop patient drawings, define the cortical/visual model, simulate and
select best-matching single-electrode phosphenes, combine them into small
sets, generate matched random phosphenes/combinations in parallel, and
compare real-vs-random statistics.

Mirrors the CURRENT (bug-fixed) state of the MATLAB Main script in this repo:
  - DO_SIMULATE_AND_COMBINE actually gates the simulate+combine block (the
    `if 0` wrapper bug found during review is already fixed upstream).
  - plot_corr_histograms / compare_real_vs_rand_stats are called once per
    drawing in the main analysis loop, not prematurely before combos exist.

vbl.enable_prescreen / vbl.prescreen_size / vbl.prescreen_thresh are set here
because the MATLAB Main script sets them, but -- confirmed by grepping
ed.m -- no function in ed.m (or ed.py) ever reads them; they are vestigial
config for an unimplemented prescreen feature. They're kept for fidelity
with the MATLAB script but have no behavioral effect.

Parallelism: MATLAB's `parfor rep = 1:vbl.n_Reps` over
ed.simulate_save_rand_drawings is ported to a ThreadPoolExecutor (when
vbl.use_thread_pool, matching MATLAB's parpool('threads')) or a
ProcessPoolExecutor (matching parpool('Processes')) over the same function.
Each rep only reads c/v/trl/tp/vbl and writes its own independent output
file, so the two pools are behaviorally equivalent here (no shared mutable
state between reps).

RNG: MATLAB's `rng(42,'twister')` is not bit-for-bit reproducible in NumPy's
Mersenne Twister (different seeding/stream algorithms), so
vbl.use_reproducible_rng only guarantees *within-Python* run-to-run
reproducibility, not cross-language numeric parity with MATLAB.
"""

from __future__ import annotations

import os
import concurrent.futures as cf

import numpy as np

import ed


def main():
    # ---------------------------
    # Setup and configuration
    # ---------------------------
    vbl = ed.setup("auto")
    vbl.enable_prescreen = 1
    vbl.prescreen_size = (64, 64)
    vbl.prescreen_thresh = 0.20

    vbl.use_thread_pool = 1 if os.name == "posix" else 0
    vbl.use_reproducible_rng = 0

    if vbl.use_reproducible_rng:
        np.random.seed(42)

    # ---------------------------
    # Load/crop patient drawings
    # ---------------------------
    # Run once if needed:
    # ed.crop_drawings(vbl)
    p_draw = ed.load_patient_drawings(vbl)

    # ---------------------------
    # Define cortical/visual model
    # ---------------------------
    c, v, trl, tp = ed.define_cortical_model(vbl)

    # ---------------------------
    # Generate and save "best" singletons
    # ---------------------------
    DO_SIMULATE_AND_COMBINE = True
    if DO_SIMULATE_AND_COMBINE:
        p_draw = ed.simulate_drawings(c, v, trl, tp, p_draw, vbl)
        for d in range(vbl.n_Drawings):
            pd = p_draw[d]
            for sd in range(pd.n_SubImages):
                target_len = pd.patient_img[sd].size
                pool_rows = pd.subimg[sd].shape[0]
                n_finite = int(np.sum(np.isfinite(pd.corr[sd])))
                print(
                    f"DRAW {vbl.dirList[d].name} d={d} sd={sd} | "
                    f"targetLen={target_len} poolRows={pool_rows} | "
                    f"finiteCorr={n_finite}/{pd.corr[sd].size}"
                )

        ed.save_simulated_drawings(p_draw, vbl)

        # ed.visualize_singletons(vbl, 0)
        ed.combine_sim_draws(vbl)

        for d in range(vbl.n_Drawings):
            ed.visualize_combinations(vbl, d)

    # ---------------------------
    # Create the random versions
    # ---------------------------
    # Remove heavy fields before passing copies to workers.
    #
    # NOTE: v.X/v.Y are intentionally NOT stripped here, unlike the MATLAB
    # Main script (which drops them via safe_rmfield(v, {'e','x','y','X','Y'})).
    # Tracing the call chain: simulate_save_rand_drawings -> generate_corticalmap
    # (recomputes c.x/c.y/c.X/c.Y and c.v.X/c.v.Y, but never touches v.X/v.Y)
    # -> generate_corticalelectricalresponse, which needs v.X's *shape* to
    # size rfmap = zeros([size(v.X), 2]) (p2p_c.m). v.X/v.Y are set once by
    # define_visualmap and never regenerated anywhere else, so stripping them
    # here would raise an undefined-field error in every worker -- this
    # looks like a latent bug in the MATLAB Main script itself (unrelated to
    # porting), not something introduced by this port. Flagged for the
    # author to confirm/fix upstream; c's heavy fields are still stripped
    # since generate_corticalmap does regenerate all of those.
    c = ed.safe_rmfield(c, ["e", "x", "y", "X", "Y", "v", "cropPix"])
    v = ed.safe_rmfield(v, ["e", "x", "y"])

    Executor = cf.ThreadPoolExecutor if vbl.use_thread_pool else cf.ProcessPoolExecutor
    max_workers = 4 if vbl.use_thread_pool else None
    with Executor(max_workers=max_workers) as ex:
        futures = [
            ex.submit(ed.simulate_save_rand_drawings, rep, c, v, trl, tp, vbl)
            for rep in range(vbl.n_Reps)
        ]
        for f in futures:
            f.result()

    print("Done generating random images")

    # ---------------------------
    # Combine random models
    # ---------------------------
    ed.combine_random_models(vbl)

    # ---------------------------
    # Analysis: real vs random
    # ---------------------------
    for d in range(vbl.n_Drawings):
        ed.show_best_combos(vbl, d)
        ed.plot_corr_histograms(vbl, d)
        ed.compare_real_vs_rand_stats(vbl, d)

    # Example for last drawing explicitly:
    last_idx = vbl.n_Drawings - 1
    ed.plot_corr_histograms(vbl, last_idx)
    ed.compare_real_vs_rand_stats(vbl, last_idx)


if __name__ == "__main__":
    main()
