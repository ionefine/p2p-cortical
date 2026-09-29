# Elche Drawings: MATLAB -> Python port notes

This directory is a faithful Python port of the Elche Drawings pipeline:
`p2p_c.m` (cortical/visual model), the `imregcode` similarity-registration
helpers, `crop_img.m`, `ed.m` (pipeline orchestration), and
`Simulate_Elche_Drawings_MainFast_5_3_2026.m` (entry-point script). The goal
throughout was **numerical fidelity to the actual MATLAB source**, not a
best-guess reimplementation -- every MATLAB source file that was needed
(including transitively) was read in full before porting it.

## File map

| Python file | Ports | Status |
|---|---|---|
| `registration.py` | `findScaleRotationNGC_if.m`, `resolveSimilarityRotationAmbiguityNGC_if.m`, `fftPadSize_if.m`, `logPolarResample_if.m`, `complexGradientImage_if.m`, `findTranslationNGC_if.m`, `normalizedGradientCorrelation_if.m`, `peakLocation2D_if.m`, `fftCorrelation2D_if.m`, `imwarp`/`affinetform2d` | Validated with a synthetic known-transform test |
| `crop_img.py` | `crop_img.m` | Validated, including the homogeneous-image edge case |
| `p2p_c.py` | `p2p_c.m` (cortical/visual model, receptive fields, temporal models) | Validated via integration smoke test (spatial path); temporal spiking path ported but not independently exercised (see below) |
| `ed.py` | `ed.m` | Validated via full-pipeline integration smoke test |
| `run_elche_drawings_main.py` | `Simulate_Elche_Drawings_MainFast_5_3_2026.m` | Validated via full-pipeline integration smoke test (both thread-pool and process-pool backends) |

Install dependencies with `pip install -r requirements.txt`.

## Running it

```
python3 run_elche_drawings_main.py
```

It will prompt "Delete all previous simulations? 1=yes, 0=no ..." (same
interactive prompt as `ed.setup` in MATLAB), then expects the working
directory (or an explicit `loc` argument to `setup`) to contain drawing
folders whose names contain "img", each with a `drawings/` subfolder holding
`<name>_edit_crop_<sd>.png` files -- identical layout to the MATLAB version.

## Key porting decisions

- **Struct semantics.** MATLAB structs are copy-on-write value types; Python
  objects are references. `p2p_c.Struct` (a thin `types.SimpleNamespace`
  subclass) is used everywhere structs were used in MATLAB. Every place the
  MATLAB code does something like `c.e = c_orig.e(eIdx)` to get a scratch
  single-electrode copy, this port uses `copy.deepcopy` on the electrode
  struct so later in-place mutation can't leak back into the shared
  original. Large read-only arrays (`X`, `Y`, `ORmap`, `cropPix`, ...) are
  intentionally shared by reference rather than deep-copied -- they're never
  mutated in place after `generate_corticalmap`, so sharing is behaviorally
  equivalent to MATLAB's copy-on-write and avoids needless duplication.

- **Indexing / array layout.** MATLAB is 1-based and column-major
  (Fortran order); NumPy defaults to 0-based, row-major (C order). Every
  index arithmetic site (e.g. `crop_img.py`'s `find(diff(...),1,'first')`
  ports, `ind2sub`-style unravelling in `registration.py`'s peak-finding)
  was translated explicitly and cross-checked against the MATLAB source line
  by line, rather than assumed to "just work" under NumPy defaults.

- **`imwarp`/`affinetform2d` conventions.** `registration.py` documents that
  `imwarp_similarity_auto`'s automatic output-bounding-box sizing is a
  principled reimplementation of MATLAB's (proprietary, undocumented)
  internal sizing logic, not guaranteed byte-identical. This matters less
  than it sounds: every call site in `ed.py` uses
  `imwarp_similarity_fixed` (a fixed/explicit `OutputView`, matching the
  patient-image target size), which has no such ambiguity and was validated
  directly.

- **File I/O.** MATLAB `.mat` files are replaced with Python `pickle`
  (`.pkl`) files, using the same directory layout (`models/`, `combos/`,
  `random_models/`, `figures/`) and filename convention
  (`vbl.fileidstr`, debug tag) as the MATLAB version. The two pipelines'
  output files are **not** directly interchangeable, but the naming *logic*
  is preserved 1:1.

- **SSIM.** Implemented as a standard 11x11-Gaussian-windowed
  structure-only SSIM (Wang et al. 2004, exponents `[0,0,1]`), matching
  MATLAB's `ssim(...,'Exponents',[0 0 1],'DynamicRange',255)` call sites in
  spirit but **not** reverse-engineered from MATLAB's internal `ssim()`
  implementation bit-for-bit. This is safe because in both pipelines the
  ssim value is only ever logged/saved, never used to select or rank
  candidates -- `fastcorr` drives every selection decision.

- **RNG.** NumPy's Mersenne Twister and MATLAB's `rng('twister')` are not
  bit-for-bit cross-compatible even with the same seed. `vbl.use_reproducible_rng`
  in `run_elche_drawings_main.py` (`np.random.seed(42)`) only guarantees
  within-Python run-to-run reproducibility, not numeric parity with a
  MATLAB run seeded the same way. Anything downstream of the random
  cortical map (`generate_corticalmap`'s `ORmap`/`ODmap`/`ONOFFmap`/`DISTmap`,
  random electrode placement, random-rep phosphene shuffles) will therefore
  differ numerically from a MATLAB run even with matched seeds -- only the
  *algorithm* is guaranteed identical, not the specific random draws.

## Preserved intentional MATLAB behaviors (confirmed with the author, not bugs)

- Random phosphenes are **not** rotated during registration
  (`simulate_save_rand_drawings` forces `r = 0`), unlike real singletons.
- Only the **first** electrode of a multi-electrode random combo is actually
  randomized (`combine_random_models` / `show_best_combos`); the rest stay
  real. Intentional: the question is "does adding electrode N+1 beat
  chance," not "is the whole combo no better than chance."
- `combine_sim_draws` / `combine_random_models` always regress against
  `patient_img[0]` (the first subimage), even when pooling candidate columns
  from every subimage.
- `norm255`'s rectify branch divides by `vmax` instead of `vmax - vmin`; a
  no-op whenever `vmin == 0`, which is guaranteed at every call site here
  (`imwarp`'s zero-fill background).

## Bugs found and fixed during the port

These were caught by integration smoke tests, root-caused against the exact
MATLAB source, and fixed (not just papered over):

1. **`p2p_c.py` / `generate_corticalelectricalresponse`.** `rfmap` was being
   written onto `c.e[idx]` instead of `v.e[idx]` (MATLAB: `v.e(idx(ii)).rfmap = ...`).
   Fixed.
2. **`ed.py` / `compare_real_vs_rand_stats`.** MATLAB's `prctile([], ...)`
   degrades gracefully to `NaN` (and MATLAB's `if [] < eps` / `if best_real_corr >= NaN`
   are simply falsy); `np.percentile([], ...)` raises `IndexError`. This is
   reachable whenever `nimg=3` has zero real or random combos (e.g.
   `vbl.n_DrawCombined < 3`, as in small/debug runs) -- the MATLAB source
   hardcodes `for nimg = 1:3` unconditionally (`ed.m:965`), so this isn't a
   debug-only corner case. Fixed by guarding the empty-array case and
   returning `NaN`, matching MATLAB's graceful degradation.
3. **`ed.py` / `generate_phosphene`.** MATLAB's `m = max(abs(img(:)))` on an
   empty `img` (a fully degenerate/uniform, e.g. no-visible-response,
   phosphene after `crop_img(img, 20)`) returns `[]`, `if [] < eps` is
   falsy (skipping the `m=1` fallback), and `img ./ []` stays empty --
   MATLAB propagates a harmless empty array rather than erroring.
   `np.max(np.abs(empty_array))` raises `ValueError`. Fixed by special-casing
   `img.size == 0` to return the empty array directly, matching MATLAB.

## Bug found in the MATLAB source (now fixed there too)

`Simulate_Elche_Drawings_MainFast_5_3_2026.m` used to strip `v.X`/`v.Y`
before dispatching the `parfor rep = 1:vbl.n_Reps` random-rep loop:

```matlab
v = ed.safe_rmfield(v, {'e','x','y','X','Y'});
...
parfor rep = 1:vbl.n_Reps
    ed.simulate_save_rand_drawings(rep, c, v, trl, tp, vbl);
end
```

Tracing the call chain: `simulate_save_rand_drawings` calls
`generate_corticalmap(c, v)`, which recomputes `c.x/c.y/c.X/c.Y` and
`c.v.X/c.v.Y` (cortical points *in visual coordinates* -- a different field
under `c.v`, not `v.X`/`v.Y`) but never touches `v.X`/`v.Y` at all
(`p2p_c.m:356-443`). `v.X`/`v.Y` (the visual field's own pixel grid) are set
exactly once, by `define_visualmap`, and never regenerated anywhere else.
`generate_corticalelectricalresponse` needs `v.X`'s *shape* directly:
`rfmap = zeros([size(v.X), 2]);` (`p2p_c.m:34`). So every parfor worker
would have hit "Reference to non-existent field 'X'" in MATLAB too -- this
was a genuine, previously-unexercised bug in the Main script, not something
introduced by porting.

**Fixed in the MATLAB source** (`Simulate_Elche_Drawings_MainFast_5_3_2026.m`):
`v = ed.safe_rmfield(v, {'e','x','y'});` -- `'X','Y'` dropped from the list,
matching `run_elche_drawings_main.py`'s already-correct behavior. `c`'s
heavy fields are still stripped (all regenerated by `generate_corticalmap`).
`run_elche_drawings_main.py` already had the correct
`safe_rmfield(v, ["e", "x", "y"])` call (validated end-to-end with both the
thread-pool and process-pool parallelism backends), so the two pipelines now
match.

## Validation performed

- `registration.py`: synthetic known-transform test (scale=1.4,
  rotation=25 deg) recovered scale=1.381, rotation=24.85 deg,
  post-alignment correlation=0.996.
- `p2p_c.py`: integration smoke test mirroring the Elche
  `define_cortical_model` + `simulate_drawings` per-electrode flow (caught
  and fixed bug #1 above).
- `crop_img.py`: unit test including the uniform-image edge case.
- `ed.py`: full smoke test through `setup`-equivalent manual construction ->
  `load_patient_drawings` (synthetic PNG test images) ->
  `define_cortical_model` (debug-scale) -> `simulate_drawings` ->
  `save_simulated_drawings` -> `combine_sim_draws` ->
  `simulate_save_rand_drawings` -> `combine_random_models` ->
  `visualize_singletons`/`visualize_combinations`/`show_best_combos`/
  `plot_corr_histograms` -> `compare_real_vs_rand_stats`; also a standalone
  `crop_drawings` test. Caught and fixed bugs #2 and #3 above.
- `run_elche_drawings_main.py`: full end-to-end smoke test (with
  `ed.setup`/`ed.define_cortical_model` monkeypatched to use tiny synthetic
  data and a trimmed electrode count for speed), exercised through both the
  `ThreadPoolExecutor` (`vbl.use_thread_pool=1`) and `ProcessPoolExecutor`
  (`vbl.use_thread_pool=0`) parallelism backends -- the latter also confirms
  `c`/`v`/`trl`/`tp`/`vbl` (all built from `p2p_c.Struct`/NumPy arrays) are
  picklable for `multiprocessing`.

## Known lower-confidence area

`vbl.trl.freq` is `NaN` throughout the Elche pipeline
(`define_cortical_model` sets `trl = Struct(freq=float("nan"))`), which
means `p2p_c.generate_phosphene` always takes the fast/spatial-only path
(`curr_trl.max_phosphene = v.e[e_idx].rfmap`) and **never** exercises the
temporal spiking model (`spike_model`, `convolve_model`, `integratefire_model`,
gamma-cascade impulse responses, etc.). Those functions were ported in full
per "port everything," and follow the MATLAB source closely, but are not
independently exercised by any test here since the Elche pipeline itself
never reaches that code path.
