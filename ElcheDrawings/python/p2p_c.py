"""
Python port of p2p_c.m — the cortical/temporal model underlying the Elche
Drawings pipeline (Fine & Boynton cortical visual prosthesis model).

Only *active* (non-commented) MATLAB methods are ported. p2p_c.m also
contains large blocks of commented-out MATLAB code (the Array_Sim_* family,
an older generate_corticalelectricalresponse/generate_corticalcell,
weibull, chronaxie, fillSymbols) which are dead code in the source and are
not translated here. Pure MATLAB-graphics functions (plotcortgrid,
plotretgrid, logx2raw, logy2raw, draw_ellipse) are also not ported, since
"porting" a MATLAB `figure`/`plot` call has no well-defined Python
equivalent without inventing a plotting API nobody asked for; everything
that produces *data* (as opposed to a figure) is ported in full.

Important: the Elche Main script sets `trl.freq = NaN` before calling
`define_trial`, which makes `generate_phosphene` take the
`isnan(trl.freq)` branch and skip the entire temporal spiking model
(spike_model -> convolve_model/integratefire_model -> slowgamma ->
nonlinearity). That temporal machinery is still ported below in full for
completeness, but it is NOT exercised by the Elche pipeline itself, so it
has had less opportunity for numerical cross-checking here than the
spatial/cortical-map code path that Elche actually depends on.

Struct convention: MATLAB structs (c, v, tp, trl, nc, and their `.e`
struct-array fields) are represented here as `types.SimpleNamespace`
objects (via the small `Struct` subclass below), so `c.pixpermm` in MATLAB
becomes `c.pixpermm` in Python too. `if ~isfield(x,'f'); x.f = val; end`
becomes `x.setdefault('f', val)`.

Randomness: no attempt is made to reproduce MATLAB's RNG stream bit-for-bit
(the algorithms differ), so stochastic parts of the model (electrode
placement, cortical map noise) will not match MATLAB run-for-run even with
"the same seed" -- only statistical/distributional equivalence should be
expected.
"""

from __future__ import annotations

import math
import types
from typing import Optional

import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import norm as _norm


# ---------------------------------------------------------------------------
# struct helper
# ---------------------------------------------------------------------------
class Struct(types.SimpleNamespace):
    """MATLAB-struct-like namespace with isfield/setdefault helpers."""

    def isfield(self, name: str) -> bool:
        return hasattr(self, name) and getattr(self, name) is not None

    def setdefault(self, name: str, value):
        if not self.isfield(name):
            setattr(self, name, value)
        return getattr(self, name)


def _field_list(struct_array, name):
    """Equivalent of MATLAB's [s.field] comma-separated-list -> array."""
    return np.array([getattr(s, name) for s in struct_array])


# ===========================================================================
# Temporal parameters / trial definition
# ===========================================================================
def define_temporalparameters(tp: Optional[Struct] = None) -> Struct:
    """Defaults based on Fine and Boynton, 2024 (Scientific Reports)."""
    if tp is None:
        tp = Struct()

    tp.setdefault("dt", 0.001 * 1e-3)
    tp.setdefault("tau1", 0.0003)
    tp.setdefault("refrac", 100)
    tp.setdefault("delta", 0.001)

    tp.setdefault("tSamp", 1000)
    tp.setdefault("tau2", 0.025)
    tp.setdefault("ncascades", 3)
    tp.setdefault("gammaflag", 1)

    tp.setdefault("adapt_on", 0)
    tp.setdefault("Na_recovery", 0.015)
    tp.setdefault("Na_strength", 500)

    tp.setdefault("spikemodel", "convolve")

    if tp.spikemodel == "convolve":
        tp.setdefault("saturation_model", "compression")
        tp.setdefault("sc_in", 0.5663)
    elif tp.spikemodel == "integratefire":
        tp.setdefault("saturation_model", "linear")
        tp.setdefault("sc_in", 1)
    else:
        raise ValueError('model variant not defined, should be "compression" or "linear"')

    if tp.saturation_model == "compression":
        tp.setdefault("power", 15.5901)
        tp.setdefault("sc_out", 10)
    elif tp.saturation_model == "sigmoid":
        tp.asymptote = 2000
        tp.e50 = 500
    elif tp.saturation_model == "normcdf":
        tp.asymptote = 1500
        tp.mean = 750
        tp.sigma = 175
    elif tp.saturation_model == "weibull":
        tp.asymptote = 1000
        tp.thresh = 600
        tp.beta = 3.5

    return tp


def generate_pt(trl: Struct, tp: Struct) -> Struct:
    if trl.isfield("on"):
        trl.lag = trl.on
        trl.dur = trl.off - trl.on

    if trl.freq is None or (isinstance(trl.freq, float) and math.isnan(trl.freq)):
        trl.pt = np.array([1.0])
    else:
        on = np.mod(trl.t, 1.0 / trl.freq) < trl.pw
        delay = trl.pw + trl.ip
        off = np.mod(trl.t - delay, 1.0 / trl.freq) < trl.pw
        tmp = trl.amp * (on.astype(float) - off.astype(float))

        lag = int(round(trl.lag / tp.dt))
        trl.pt = np.zeros(lag + len(tmp))
        trl.pt[lag : lag + len(tmp)] = tmp

        trl.t = np.arange(0, trl.dur + trl.lag, tp.dt)
        trl.t = trl.t[:-1]

    if trl.dur < trl.simdur:
        n_needed = int(round(trl.simdur / tp.dt))
        if n_needed > len(trl.pt):
            trl.pt = np.concatenate([trl.pt, np.zeros(n_needed - len(trl.pt))])
        trl.t = np.arange(0, trl.simdur, tp.dt)

    return trl


def define_trial(tp: Struct, trl: Optional[Struct] = None) -> Struct:
    if trl is None:
        trl = Struct()

    trl.setdefault("e", 1)
    trl.setdefault("dur", 1000 * 1e-3)
    trl.setdefault("simdur", 3)

    trl.t = np.arange(0, trl.dur - tp.dt, tp.dt)
    # np.arange with a stop that's an exact multiple of the step can behave
    # differently from MATLAB's colon operator at the floating-point edge;
    # match MATLAB's `0:dt:dur-dt` element count exactly:
    n = int(round((trl.dur - tp.dt) / tp.dt)) + 1
    trl.t = np.arange(n) * tp.dt

    trl.setdefault("pw", 0.1 * 1e-3)
    trl.setdefault("ip", trl.pw)
    trl.setdefault("lag", 2 * trl.pw)
    trl.setdefault("order", 1)
    trl.setdefault("freq", 60)
    trl.setdefault("amp", 100)

    trl = generate_pt(trl, tp)

    trl.CperTrial = (trl.amp / 1000) * trl.dur * trl.freq * trl.pw * 1e3
    trl.CperPulse = trl.pw * trl.amp / 1000

    return trl


# ===========================================================================
# Cortex / visual map definitions
# ===========================================================================
def define_cortex(c: Optional[Struct] = None) -> Struct:
    if c is None:
        c = Struct()

    c.setdefault("efthr", 0.05)
    c.setdefault("animal", "human")

    if c.animal == "human":
        c.k = 15
        c.setdefault("a", 0.5)
        c.shift = c.k * math.log(c.a)
        c.setdefault("squish", 1)
        c.setdefault("cortexHeight", [-40, 40])
        c.setdefault("cortexLength", [-5, 80])
        c.setdefault("pixpermm", 8)
    elif c.animal == "macaque":
        c.k = 5
        c.squish = 1
        c.a = 0.3
        c.shift = c.k * math.log(c.a)
        c.setdefault("cortexHeight", [-20, 20])
        c.setdefault("cortexLength", [-5, 30])
        c.setdefault("pixpermm", 8)
    elif c.animal == "mouse":
        raise NotImplementedError("Sorry no model for mouse yet")
    else:
        raise ValueError(f"Unrecognized c.animal: {c.animal}")

    c.setdefault("rfmodel", "ringach")
    c.setdefault("rfsizemodel", "keliris")

    if c.animal == "human":
        if c.rfsizemodel == "keliris":
            c.slope = 0.08
            c.intercept = 0.16
            c.min = 0
        elif c.rfsizemodel == "bosking":
            c.slope = 0.2620 / 2
            c.intercept = 0.1787 / 2
            c.min = 0
        elif c.rfsizemodel == "winawer":
            c.slope = 0.1667
            c.min = 1.11
            c.intercept = 0.0721
    elif c.animal == "macaque":
        if c.rfsizemodel == "keliris":
            c.slope = 0.08
            c.intercept = 0.16
            c.min = 0
        else:
            c.slope = 0.06
            c.intercept = 0.42
            c.min = 0

    c.setdefault("ar", 0.25)
    c.setdefault("delta", 2)
    c.setdefault("onoff_ratio", 0.8)
    c.setdefault("sig", 0.5)

    if c.animal == "human":
        c.ODsize = 0.863
        c.filtSz = 3
    elif c.animal == "macaque":
        c.ODsize = 0.531
        c.filtSz = 1.85
    elif c.animal == "mouse":
        c.ODsize = np.nan
        c.filtSz = np.nan

    c.gridColor = [1, 1, 0]
    return c


def define_visualmap(v: Optional[Struct] = None) -> Struct:
    if v is None:
        v = Struct()

    v.setdefault("visfieldHeight", [-30, 30])
    v.setdefault("visfieldWidth", [-30, 30])
    v.setdefault("pixperdeg", 7)
    v.setdefault("drawthr", 1)

    nx = int((v.visfieldWidth[1] - v.visfieldWidth[0]) * v.pixperdeg)
    ny = int((v.visfieldHeight[1] - v.visfieldHeight[0]) * v.pixperdeg)
    v.x = np.linspace(v.visfieldWidth[0], v.visfieldWidth[1], nx)
    v.y = np.linspace(v.visfieldHeight[0], v.visfieldHeight[1], ny)

    v.X, v.Y = np.meshgrid(v.x, v.y)

    v.setdefault("angList", list(range(-90, 91, 45)))
    v.setdefault("eccList", [1, 2, 3, 5, 8, 13, 21, 34])
    v.gridColor = [1, 1, 0]
    v.n = 201

    return v


# ===========================================================================
# Coordinate transforms (Schwartz log-polar cortical map)
# ===========================================================================
def v2c_cplx(c: Struct, z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=np.complex128).copy()
    lvf = np.real(z) < 0
    z[lvf] = z[lvf] * np.exp(-1j * np.pi)
    w = (c.k * np.log(z + c.a)) - c.shift
    w[~lvf] = w[~lvf] * np.exp(-1j * np.pi)
    w = np.real(w) + c.squish * 1j * np.imag(w)
    return w


def v2c_real(c: Struct, vx, vy):
    z = v2c_cplx(c, np.asarray(vx) + 1j * np.asarray(vy))
    return np.real(z), np.imag(z)


def isValidCortex(c: Struct, cx, cy):
    cx = np.asarray(cx, dtype=float)
    cy = np.asarray(cy, dtype=float)
    vyb = c.a * np.tan(np.abs(cy) / c.k)
    cxb = c.k * np.log(np.sqrt(vyb**2 + c.a**2)) - c.shift
    return (np.abs(cx) > cxb) & (np.abs(cy) < c.k * np.pi / 2)


def c2v_cplx(c: Struct, z: np.ndarray):
    z = np.asarray(z, dtype=np.complex128).copy()
    z = np.real(z) + 1j * np.imag(z) / c.squish

    lvf = np.real(z) > 0
    z[~lvf] = z[~lvf] * np.exp(-1j * np.pi)
    w = np.exp((z + c.shift) / c.k) - c.a
    w[lvf] = w[lvf] * np.exp(-1j * np.pi)

    ok = isValidCortex(c, np.real(z), np.imag(z))
    w = np.where(ok, w, np.nan + 1j * np.nan)
    return w, ok


def c2v_real(c: Struct, cx, cy):
    z, ok = c2v_cplx(c, np.asarray(cx) + 1j * np.asarray(cy))
    return np.real(z), np.imag(z), ok


# ===========================================================================
# Electrodes
# ===========================================================================
def define_electrodes(c: Struct, v: Struct):
    """Electrodes given in visual-field coordinates -> placed on cortex."""
    for e in v.e:
        if not hasattr(e, "ang"):
            ang, ecc = _cart2pol(e.x, e.y)
            e.ang = ang * 180 / np.pi
            e.ecc = ecc

    for e in c.e:
        if not e.isfield("radius") if isinstance(e, Struct) else not hasattr(e, "radius"):
            e.radius = 500 / 1000
        if not hasattr(e, "shape"):
            e.shape = "round"

    ang_all = np.array([e.ang for e in v.e]) * np.pi / 180
    ecc_all = np.array([e.ecc for e in v.e])
    x_all, y_all = _pol2cart(ang_all, ecc_all)
    for i, e in enumerate(v.e):
        e.x = x_all[i]
        e.y = y_all[i]

    radii = np.array([e.radius for e in c.e])
    areas = np.pi * radii**2
    vx = np.array([e.x for e in v.e])
    vy = np.array([e.y for e in v.e])
    cx, cy = v2c_real(c, vx, vy)
    for i, e in enumerate(c.e):
        e.area = areas[i]
        e.x = cx[i]
        e.y = cy[i]

    return c, v


def c2v_define_electrodes(c: Struct, v: Struct) -> Struct:
    """Electrodes given on cortex -> projected into visual-field coordinates."""
    cx = np.array([e.x for e in c.e])
    cy = np.array([e.y for e in c.e])
    vx, vy, _ok = c2v_real(c, cx, cy)
    angs, eccs = _cart2pol(vx, vy)
    angs = angs * 180 / np.pi
    for i, e in enumerate(v.e):
        e.x = vx[i]
        e.y = vy[i]
        e.ang = angs[i]
        e.ecc = eccs[i]
    return v


def _cart2pol(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.arctan2(y, x), np.hypot(x, y)


def _pol2cart(theta, r):
    theta = np.asarray(theta, dtype=float)
    r = np.asarray(r, dtype=float)
    return r * np.cos(theta), r * np.sin(theta)


# ===========================================================================
# Electric field
# ===========================================================================
def generate_ef(c: Struct):
    c.setdefault("emodel", "Tehovnik")

    for e in c.e:
        R = np.sqrt((c.X - e.x) ** 2 + (c.Y - e.y) ** 2)
        Rd = R - e.radius
        Rd = np.where(Rd < 0, 0, Rd)

        if c.emodel == "Tehovnik":
            c.setdefault("I_0", 1)
            c.setdefault("I_k", 6.75)
            pt_ef = c.I_0 / (1 + c.I_k * Rd**2)
        else:
            raise ValueError("electric field model not specified in code")

        scaled = 255.0 * pt_ef / np.max(pt_ef)
        e.ef = np.clip(np.round(scaled), 0, 255).astype(np.uint8)

    return c


# ===========================================================================
# Cortical map generation (orientation, ocular dominance, on/off, RF size)
# ===========================================================================
def generate_corticalmap(c: Struct, v: Struct):
    nx = int((max(c.cortexLength) - min(c.cortexLength)) * c.pixpermm)
    ny = int((max(c.cortexHeight) - min(c.cortexHeight)) * c.pixpermm)
    c.x = np.linspace(min(c.cortexLength), max(c.cortexLength), nx)
    c.y = np.linspace(min(c.cortexHeight), max(c.cortexHeight), ny)
    c.X, c.Y = np.meshgrid(c.x, c.y)
    sz = c.X.shape

    # --- orientation & ocular dominance maps: bandpass-filtered white noise
    # (Rojer & Schwartz, 1990) ---
    Z = np.exp(1j * np.random.rand(*sz) * 2 * np.pi)

    freq = 1.0 / c.ODsize
    filtPix = int(np.ceil(c.filtSz * c.pixpermm))
    fx = np.linspace(-c.filtSz / 2, c.filtSz / 2, filtPix)
    Xf, Yf = np.meshgrid(fx, fx)
    Rf = np.sqrt(Xf**2 + Yf**2)
    FILT = np.exp(-(Rf**2) / c.sig**2) * np.cos(2 * np.pi * freq * Rf)

    W = fftconvolve(Z, FILT, mode="same")
    c.ORmap = np.angle(W)
    WX = np.gradient(W, axis=1)
    Gx = np.angle(WX)
    c.ODmap = _norm.cdf(Gx * c.sig)

    # --- on/off & distance maps ---
    freq = (1.0 / c.ODsize) * 2
    filtPix = int(np.ceil(c.filtSz / 2 * c.pixpermm))
    fx = np.linspace(-c.filtSz / 2, c.filtSz / 2, filtPix)
    Xf, Yf = np.meshgrid(fx, fx)
    Rf = np.sqrt(Xf**2 + Yf**2)
    FILT = np.exp(-(Rf**2) / c.sig**2) * np.cos(2 * np.pi * freq * Rf)

    W = fftconvolve(Z, FILT, mode="same")
    u = np.angle(W) / np.pi
    tmp = np.zeros_like(u)
    pos = u > 0
    neg = u < 0
    tmp[pos] = -np.log(u[pos]) / c.delta
    tmp[neg] = np.log(-u[neg]) / c.delta
    c.DISTmap = tmp
    WX = np.gradient(W, axis=1)
    Gx = np.angle(WX)
    c.ONOFFmap = _norm.cdf(Gx * c.sig)

    # --- angle / eccentricity maps ---
    vX, vY, _ok = c2v_real(c, c.X, c.Y)
    c.v = Struct()
    c.v.X = vX
    c.v.Y = vY
    ANG, ECC = _cart2pol(vX, vY)
    c.v.ANG = ANG * 180 / np.pi
    c.v.ECC = ECC

    # grid lines (for plotting only in MATLAB; kept as data here)
    rho = np.linspace(0, max(v.eccList), v.n)
    ang = np.array(v.angList) * np.pi / 180
    v.zAng = rho[:, None] * np.exp(1j * ang[None, :])
    c.v.gridAng = v2c_cplx(c, v.zAng)

    eccList = np.array(v.eccList, dtype=float)
    ang2 = np.linspace(-90, 90, v.n) * np.pi / 180
    v.zEcc = (eccList[:, None] * np.exp(1j * ang2[None, :])).T
    c.v.gridEcc = v2c_cplx(c, v.zEcc)

    c.RFsizemap = np.maximum(c.slope * np.abs(c.v.ECC) + c.intercept, c.min)

    c.cropPix = c.v.ANG.copy()
    ecc_limit = max(max(v.visfieldHeight), max(v.eccList))
    c.cropPix[c.v.ECC > ecc_limit] = np.nan
    if abs(min(c.cortexLength)) < abs(max(c.cortexLength)):
        c.cropPix[c.X < 0] = np.nan
    elif abs(min(c.cortexLength)) > abs(max(c.cortexLength)):
        c.cropPix[c.X > 0] = np.nan

    return c, v


# ===========================================================================
# Receptive fields / cortical electrical response
# ===========================================================================
def _flat_maps_cache(c: Struct):
    """Column-major-flattened views of the per-pixel maps generate_corticalcell
    indexes into. These maps (c.v.X/Y, c.ODmap, c.ORmap, c.RFsizemap,
    c.DISTmap, c.ONOFFmap) don't change across the pixel loop in
    generate_corticalelectricalresponse, so flattening them once and reusing
    the result (instead of re-flattening -- i.e. reallocating a full copy --
    on every single pixel call) is a pure performance fix with no numerical
    effect: same values, just computed once instead of thousands of times."""
    return {
        "vX": c.v.X.flatten(order="F"),
        "vY": c.v.Y.flatten(order="F"),
        "OD": c.ODmap.flatten(order="F"),
        "OR": c.ORmap.flatten(order="F"),
        "RFsize": c.RFsizemap.flatten(order="F"),
        "DIST": c.DISTmap.flatten(order="F"),
        "ONOFF": c.ONOFFmap.flatten(order="F"),
    }


def generate_corticalcell(ef, pix_flat_index: int, c: Struct, v: Struct, flat=None):
    """RF for a single cortical pixel (flattened MATLAB column-major index).

    Returns an array shaped (ny, nx, 2) for 'scoreboard'/'smirnakis', or
    (ny, nx, 2, 2) for 'ringach' (dims: eye, on/off).

    `flat`, if given, is a cache from `_flat_maps_cache(c)` -- avoids
    re-flattening c's per-pixel maps on every call (see that function's
    docstring). Falls back to flattening internally when not provided, so
    the standalone call signature is unchanged.
    """
    if flat is None:
        flat = _flat_maps_cache(c)

    x0 = flat["vX"][pix_flat_index]
    y0 = flat["vY"][pix_flat_index]
    od = flat["OD"][pix_flat_index]
    theta = np.pi - flat["OR"][pix_flat_index]
    sigma_x = flat["RFsize"][pix_flat_index] * c.ar
    sigma_y = flat["RFsize"][pix_flat_index]

    if c.rfmodel == "scoreboard":
        G = ef * np.exp(-((v.X - x0) ** 2 / 0.0001 + (v.Y - y0) ** 2 / 0.00001))
        RF = np.stack([G, G], axis=-1)
        return RF

    aa = np.cos(theta) ** 2 / (2 * sigma_x**2) + np.sin(theta) ** 2 / (2 * sigma_y**2)
    bb = -np.sin(2 * theta) / (4 * sigma_x**2) + np.sin(2 * theta) / (4 * sigma_y**2)
    cc = np.sin(theta) ** 2 / (2 * sigma_x**2) + np.cos(theta) ** 2 / (2 * sigma_y**2)

    if c.rfmodel == "smirnakis":
        G = ef * np.exp(-(aa * (v.X - x0) ** 2 + 2 * bb * (v.X - x0) * (v.Y - y0) + cc * (v.Y - y0) ** 2))
        G = G / np.sum(G)
        RF = np.stack([od * G, (1 - od) * G], axis=-1)
        return RF

    if c.rfmodel == "ringach":
        tmp = np.exp(-(aa * (v.X - x0) ** 2 + 2 * bb * (v.X - x0) * (v.Y - y0) + cc * (v.Y - y0) ** 2))
        A = np.sqrt(np.sum(tmp > 0.2) / v.pixperdeg**2)
        pixNum_1based = pix_flat_index + 1  # DISTmap column-major index below
        d = flat["DIST"][pix_flat_index] * A

        x_off = x0 + (d / 2) * np.cos(theta)
        y_off = y0 - (d / 2) * np.sin(theta)
        x_on = x0 - (d / 2) * np.cos(theta)
        y_on = y0 + (d / 2) * np.sin(theta)

        wplus = flat["ONOFF"][pix_flat_index]
        wminus = 1 - wplus

        hplus_on = np.exp(-(aa * (v.X - x_on) ** 2 + 2 * bb * (v.X - x_on) * (v.Y - y_on) + cc * (v.Y - y_on) ** 2))
        hplus_on = hplus_on / np.max(np.abs(hplus_on))
        hplus_on = hplus_on * wplus

        hplus_off = np.exp(-(aa * (v.X - x_off) ** 2 + 2 * bb * (v.X - x_off) * (v.Y - y_off) + cc * (v.Y - y_off) ** 2))
        hplus_off = hplus_off / np.abs(np.max(hplus_off))
        hplus_off = hplus_off * wminus * c.onoff_ratio

        hminus_off = -0.4 * hplus_off
        hminus_on = -0.4 * hplus_on

        excitatory = hplus_on - c.onoff_ratio * hplus_off
        inhibitory = hminus_on - c.onoff_ratio * hminus_off

        RF = np.zeros(v.X.shape + (2, 2))
        RF[:, :, 0, 0] = od * ef * excitatory
        RF[:, :, 1, 0] = (1 - od) * ef * excitatory
        RF[:, :, 0, 1] = od * ef * inhibitory
        RF[:, :, 1, 1] = (1 - od) * inhibitory
        return RF

    raise ValueError("c.rfmodel model not recognized")


def generate_corticalelectricalresponse(c: Struct, v: Struct):
    """Sum of weighted receptive fields activated by each electrode,
    normalized so the max (absolute value) is 1.

    Note: electrode position/electric-field data is read from c.e, but the
    resulting rfmap is written onto the corresponding v.e entry -- this
    matches MATLAB's `v.e(idx(ii)).rfmap = ...` (c and v are parallel
    per-electrode arrays of the same length, in the same order).
    """
    c.setdefault("rfmodel", "ringach")

    cropFlat = c.cropPix.flatten(order="F")
    flat = _flat_maps_cache(c)

    if len(v.e) != len(c.e):
        raise ValueError("c.e and v.e must be parallel arrays of the same length")

    for idx, ce in enumerate(c.e):
        if (
            np.min(ce.x) - 1 < np.min(c.X)
            or np.max(ce.x) + 1 > np.max(c.X)
            or np.min(ce.y) - 1 < np.min(c.Y)
            or np.max(ce.y) + 1 > np.max(c.Y)
        ):
            raise ValueError("Electrode is either outside or too close to the edge of the cortical sheet")

        rfmap = np.zeros(v.X.shape + (2,))
        ef_flat = ce.ef.flatten(order="F").astype(np.float64)

        valid_mask = ~np.isnan(cropFlat) & (np.abs(ef_flat) > c.efthr * 255)
        valid_pix = np.flatnonzero(valid_mask)
        ct = len(valid_pix)

        for p in valid_pix:
            RF = generate_corticalcell(ef_flat[p], p, c, v, flat=flat)
            if RF.ndim == 4:
                RF = RF[:, :, :, 0]  # ignore inhibitory component, matches MATLAB squeeze(RF(:,:,:,1))
            rfmap[:, :, 0] += RF[:, :, 0]
            rfmap[:, :, 1] += RF[:, :, 1]

        if ct < ce.radius * 10:
            if ct == 0:
                print("WARNING! No pixels passed ef threshold.")
            else:
                print("WARNING! Very few pixels passed ef threshold.")
            print(" try checking the following:")
            print("lowering c.efthr or increase stimulation intensity")
            print("checking location of electrodes relative visual map")
            print("check the sampling resolution of cortex is not too low")

        v.e[idx].rfmap = rfmap / np.max(np.abs(rfmap))

    return v, c


# ===========================================================================
# Ellipse fitting
# ===========================================================================
def fit_ellipse_to_phosphene(img: np.ndarray, v: Struct) -> Struct:
    """Image-moment ellipse fit (mirrors p2p_c.fit_ellipse_to_phosphene)."""
    img = np.asarray(img, dtype=np.float64)
    M00 = np.sum(img)
    M10 = np.sum(v.X * img)
    M01 = np.sum(v.Y * img)
    M11 = np.sum(v.X * v.Y * img)
    M20 = np.sum(v.X**2 * img)
    M02 = np.sum(v.Y**2 * img)

    p = Struct()
    with np.errstate(invalid="ignore", divide="ignore"):
        p.x0 = M10 / M00
        p.y0 = M01 / M00
        mu20 = M20 / M00 - p.x0**2
        mu02 = M02 / M00 - p.y0**2
        mu11 = M11 / M00 - p.x0 * p.y0
    a = (mu20 + mu02) / 2
    b = 0.5 * np.sqrt(4 * mu11**2 + (mu20 - mu02) ** 2)
    lambda_1 = a + b
    lambda_2 = a - b
    p.theta = -0.5 * np.arctan2(2 * mu11, mu20 - mu02)
    p.sigma_x = 2 * np.sqrt(lambda_1)
    p.sigma_y = 2 * np.sqrt(lambda_2)
    return p


def Gauss_2D(v: Struct, x0, y0, theta, sigma_x, sigma_y):
    aa = np.cos(theta) ** 2 / (2 * sigma_x**2) + np.sin(theta) ** 2 / (2 * sigma_y**2)
    bb = -np.sin(2 * theta) / (4 * sigma_x**2) + np.sin(2 * theta) / (4 * sigma_y**2)
    cc = np.sin(theta) ** 2 / (2 * sigma_x**2) + np.cos(theta) ** 2 / (2 * sigma_y**2)
    return np.exp(-(aa * (v.X - x0) ** 2 + 2 * bb * (v.X - x0) * (v.Y - y0) + cc * (v.Y - y0) ** 2))


# ===========================================================================
# Phosphene generation
# ===========================================================================
def generate_phosphene(v: Struct, tp: Struct, trl: Struct):
    """Generate a phosphene for every electrode in v.e using a single trial
    definition `trl`. Returns a list of trial structs (one per electrode).

    Note: for the Elche pipeline, trl.freq is NaN, so this always takes the
    fast path (curr_trl.max_phosphene = v.e[e].rfmap) and never touches the
    temporal spiking model below.
    """
    n_electrodes = len(v.e)
    trl_array = []

    freq_is_nan = trl.freq is None or (isinstance(trl.freq, float) and math.isnan(trl.freq))

    for e_idx in range(n_electrodes):
        curr_trl = Struct(**vars(trl))
        curr_trl.e = e_idx

        if freq_is_nan:
            curr_trl.max_phosphene = v.e[e_idx].rfmap
            curr_trl.spikeStrength = 1
        else:
            tmp = tp.saturation_model
            tp.saturation_model = "linear"
            curr_trl = spike_model(tp, curr_trl)
            curr_trl.max_phosphene = v.e[e_idx].rfmap * curr_trl.max_temporal_response
            tp.saturation_model = tmp
            curr_trl.max_phosphene = nonlinearity(tp, curr_trl.max_phosphene)

        curr_trl.ellipse = []
        if curr_trl.max_phosphene is not None and np.size(curr_trl.max_phosphene):
            for eye in range(2):
                p = fit_ellipse_to_phosphene(curr_trl.max_phosphene[:, :, eye] > v.drawthr, v)
                curr_trl.ellipse.append(p)
        else:
            curr_trl.sim_brightness = None

        trl_array.append(curr_trl)

    return trl_array


def generate_phosphene2(v: Struct, tp: Struct, trl: Struct):
    """Single-electrode phosphene at trl.e, finding the brightest moment."""
    freq_is_nan = trl.freq is None or (isinstance(trl.freq, float) and math.isnan(trl.freq))
    if freq_is_nan:
        trl.max_phosphene = v.e[trl.e].rfmap
        trl.spikeStrength = 1
    else:
        trl = convolve_model(tp, trl)
        trl.max_phosphene = v.e[trl.e].rfmap * np.max(trl.spikeStrength)

    trl.sim_area = (1 / v.pixperdeg**2) * np.sum(trl.max_phosphene > v.drawthr) / 2
    trl.ellipse = []
    if np.size(trl.max_phosphene):
        for eye in range(2):
            p = fit_ellipse_to_phosphene(trl.max_phosphene[:, :, eye] > v.drawthr, v)
            trl.ellipse.append(p)
        beta = 6
        trl.sim_brightness = (1 / v.pixperdeg**2) * np.sum(trl.max_phosphene**beta) ** (1 / beta)
    else:
        trl.sim_brightness = None
    return trl, v


def generate_phosphene_multiple(v: Struct, tp: Struct, trl_array):
    for i, trl in enumerate(trl_array):
        freq_is_nan = trl.freq is None or (isinstance(trl.freq, float) and math.isnan(trl.freq))
        if freq_is_nan:
            trl.max_phosphene = v.e[trl.e].rfmap
            trl.spikeStrength = 1
        else:
            trl = convolve_model(tp, trl)
            trl.max_phosphene = v.e[trl.e].rfmap * np.max(trl.spikeStrength)

        trl.sim_area = (1 / v.pixperdeg**2) * np.sum(trl.max_phosphene > v.drawthr) / 2
        if np.size(trl.max_phosphene):
            trl.ellipse = []
            for eye in range(2):
                p = fit_ellipse_to_phosphene(trl.max_phosphene[:, :, eye] > v.drawthr, v)
                trl.ellipse.append(p)
            beta = 6
            trl.sim_brightness = (1 / v.pixperdeg**2) * np.sum(trl.max_phosphene**beta) ** (1 / beta)
        else:
            trl.sim_brightness = None
        trl_array[i] = trl
    return trl_array, v


# ===========================================================================
# Temporal spiking model (NOT exercised by the Elche pipeline; trl.freq=NaN
# there bypasses all of this. Ported in full for completeness.)
# ===========================================================================
def gamma(n, k, t):
    t = np.asarray(t, dtype=np.float64)
    y = (t / k) ** (n - 1) * np.exp(-t / k) / (k * math.factorial(int(n) - 1))
    y = np.where(t < 0, 0, y)
    return y


def nonlinearity(tp: Struct, x):
    x = np.asarray(x, dtype=np.float64)
    sc_in = tp.get("sc_in", 1) if hasattr(tp, "get") else getattr(tp, "sc_in", 1)

    model = tp.saturation_model
    if model == "sigmoid":
        return sc_in * x**tp.power / (x**tp.power + tp.sigma**2)
    elif model == "normcdf":
        y = _norm.cdf(x, loc=tp.mean, scale=tp.sigma)
        return np.where(y < 0, 0, y)
    elif model == "weibull":
        raise NotImplementedError("weibull nonlinearity: source weibull() is commented out in p2p_c.m")
    elif model == "power":
        return sc_in * x**tp.power
    elif model == "exp":
        return sc_in * x**tp.k
    elif model == "compression":
        return tp.sc_in * (tp.power * np.tanh((x * tp.sc_in) / tp.power))
    elif model == "linear":
        return tp.sc_in * x
    else:
        raise ValueError(f"Unrecognized tp.saturation_model: {model}")


def slowgamma(tp: Struct, trl: Struct) -> Struct:
    h = gamma(tp.ncascades, tp.tau2, trl.t)
    csum = np.cumsum(h) * tp.dt
    tid = np.flatnonzero(csum >= 0.998)
    if tid.size:
        h = h[: tid[0] + 1]
        trl.imp_resp = h
    else:
        print("Warning: gamma hdr might not have a long enough time vector")

    impFrames = np.arange(len(trl.imp_resp))
    resp = np.zeros(len(trl.t) + len(trl.imp_resp))

    if tp.spikemodel == "convolve":
        interspike = np.diff(trl.t[trl.spikeWhen], prepend=trl.t[trl.spikeWhen][0] if len(trl.spikeWhen) else 0)
        if len(trl.spikeWhen):
            interspike[0] = trl.t[trl.spikeWhen[0]]
        trl.spikes_norefrac = trl.spikeStrength
        trl.spikeStrength = trl.spikeStrength * (1 - np.exp(-tp.refrac * (interspike + tp.delta)))

    for i in range(len(trl.spikeStrength)):
        idx = np.flatnonzero(trl.t > trl.t[trl.spikeWhen[i]])
        if idx.size == 0:
            continue
        id0 = idx[0]
        resp[id0 + impFrames] += trl.spikeStrength[i] * trl.imp_resp

    trl.resp = resp
    return trl


def convolve_model(tp: Struct, trl: Struct) -> Struct:
    t = trl.t[:: tp.tSamp] if (tp.isfield("tSamp") and tp.tSamp != 1) else trl.t

    trl.spikeWhen = np.flatnonzero(np.diff(trl.pt)) + 1

    n_events = len(trl.spikeWhen)
    Rtmp = np.zeros(n_events)
    wasRising = False
    spikeId = np.zeros(n_events, dtype=bool)

    for i in range(n_events - 1):
        delta = trl.t[trl.spikeWhen[i + 1]] - trl.t[trl.spikeWhen[i]]
        Rtmp[i + 1] = trl.pt[trl.spikeWhen[i]] * tp.tau1 * (1 - np.exp(-delta / tp.tau1)) + Rtmp[i] * np.exp(
            -delta / tp.tau1
        )
        if Rtmp[i + 1] < Rtmp[i] and wasRising:
            spikeId[i] = True
            wasRising = False
        else:
            wasRising = True

    R1 = Rtmp * 1000
    R1 = np.where(R1 < 0, 0, R1)

    if R1.size == 0:
        trl.resp = np.zeros(len(t))
        trl.resp_lin = trl.resp.copy()
        trl.tt = t
        trl.max_temporal_response = 0
        trl.spikeWhen = np.array([np.nan])
        trl.imp_resp = np.array([np.nan])
        trl.spikeStrength = np.array([np.nan])
        trl.spikes_norefrac = np.array([np.nan])
    else:
        spikeId = spikeId & (R1 > 0)
        trl.spikeStrength = R1[spikeId]
        trl.spikeWhen = trl.spikeWhen[spikeId]

    return trl


def define_integratefirecells(arg=None):
    defaults = dict(gl=10, El=-58, vt=-50, delT=2, a=2, tauW=0.15, b=100, vreset=-46, vpeak=0, tARP=0)

    if isinstance(arg, list) and arg and isinstance(arg[0], Struct):
        nc = arg
        for cell in nc:
            for k, dv in defaults.items():
                cell.setdefault(k, dv)
            if not cell.isfield("R"):
                cell.R = 0.25 - math.log(np.random.rand()) / 2
        return nc

    nCells = arg if isinstance(arg, int) else 1
    nc = []
    for _ in range(nCells):
        cell = Struct(**defaults)
        cell.R = 0.25 - math.log(np.random.rand()) / 2
        nc.append(cell)
    return nc


def Na_adapt(tp: Struct, trl: Struct) -> Struct:
    h = gamma(tp.ncascades, tp.Na_recovery, trl.t)
    csum = np.cumsum(h) * tp.dt
    tid = np.flatnonzero(csum >= 0.999)
    if tid.size:
        h = h[: tid[0] + 1]
        trl.imp_resp = h
    else:
        print("Warning: gamma hdr might not have a long enough time vector")
    trl.adapt = tp.Na_strength * np.convolve(np.maximum(trl.pt, 0), h, mode="same")
    return trl


def integratefire_model(tp: Struct, trl: Struct, nc):
    nt = len(trl.t)
    dt = tp.dt
    Iinput = trl.pt * 1000

    trl.indexFire = np.full((len(trl.pt), len(nc)), np.nan)

    for cc, cell in enumerate(nc):
        nARP = int(round(cell.tARP / dt))

        d = 0
        v_ = np.zeros(nt)
        v_[0] = cell.El
        w = np.zeros(nt)
        tmp = np.zeros(nt)
        Na_adapt_ = np.ones(nt)

        for timei in range(nt - 1):
            v_[timei + 1] = v_[timei] + dt * (
                -cell.gl * (v_[timei] - cell.El)
                + cell.gl * cell.delT * np.exp((v_[timei] - cell.vt) / cell.delT)
                + (Na_adapt_[timei] * Iinput[timei])
                - w[timei]
            ) * cell.R
            if d > 0:
                v_[timei + 1] = cell.vreset

            w[timei + 1] = w[timei] + dt * (cell.a * (v_[timei] - cell.El) - w[timei]) / cell.tauW
            tmp[timei + 1] = tmp[timei] + dt * max(Iinput[timei], 0) - (dt / tp.Na_recovery) * tmp[timei]
            Na_adapt_[timei + 1] = 1 / (1 + tp.Na_strength * tmp[timei + 1])

            if v_[timei + 1] >= cell.vpeak:
                v_[timei] = cell.vpeak
                v_[timei + 1] = cell.vreset
                w[timei + 1] = w[timei + 1] + cell.b
                d = nARP
            d = max(d - 1, 0)

        trl.indexFire[:, cc] = v_ == cell.vpeak

    if len(nc) > 1:
        fired_any = np.sum(trl.indexFire, axis=1)
        trl.spikeWhen = np.flatnonzero(fired_any)
        trl.spikeStrength = np.sum(trl.indexFire[trl.spikeWhen, :], axis=1) / len(nc)
    else:
        trl.spikeWhen = np.flatnonzero(trl.indexFire[:, 0])
        trl.spikeStrength = trl.indexFire[trl.spikeWhen, 0]

    trl.Na_adapt = Na_adapt_
    return trl, nc, None


def spike_model(tp: Struct, trl: Struct, nc=None):
    tp.setdefault("spikemodel", "convolve")

    if tp.spikemodel == "convolve":
        trl = convolve_model(tp, trl)
    elif tp.spikemodel == "integratefire":
        if nc is None:
            nc = define_integratefirecells(1)
        elif isinstance(nc, int):
            nc = define_integratefirecells(nc)
        trl, nc, _ = integratefire_model(tp, trl, nc)

    if tp.gammaflag:
        trl = slowgamma(tp, trl)
    else:
        trl.resp = np.zeros(len(trl.pt))
        trl.resp[trl.spikeWhen] = trl.spikeStrength

    trl.resp = nonlinearity(tp, trl.resp)
    trl.max_temporal_response = np.max(trl.resp)
    trl.resp = trl.resp[: len(trl.pt)]
    if np.max(trl.resp) > trl.max_temporal_response:
        raise ValueError("trial has been cut off before peak response. Increase trl.simdur")

    return trl if nc is None else (trl, nc)


def find_threshold(trl: Struct, tp: Struct, nc_or_ncells=None, verbose: int = 0):
    if nc_or_ncells is None:
        nc = define_integratefirecells(1)
    elif isinstance(nc_or_ncells, list):
        nc = nc_or_ncells
    elif tp.spikemodel == "integratefire":
        nc = define_integratefirecells(nc_or_ncells)
    else:
        nc = None

    tp.setdefault("nReps", 12)

    resp = 0
    if verbose:
        print("finding upper amplitude for search")
    trl.amp = 0.5

    while resp < tp.thresh_resp:
        trl.amp = trl.amp * 2
        trl = define_trial(tp, trl)
        if tp.spikemodel == "integratefire":
            trl = spike_model(tp, trl, nc)
        elif tp.spikemodel == "convolve":
            trl = spike_model(tp, trl)
        else:
            raise ValueError("spiking model not recognized")
        resp = np.max(trl.resp)
        print(f"amp = {trl.amp} resp = {resp} thresh = {tp.thresh_resp}")

    hi = trl.amp
    lo = trl.amp / 2 if trl.amp > 1 else 0

    for _ in range(tp.nReps):
        trl.amp = (hi + lo) / 2
        trl = define_trial(tp, trl)
        if tp.spikemodel == "integratefire":
            trl = spike_model(tp, trl, nc)
        elif tp.spikemodel == "convolve":
            trl = spike_model(tp, trl)
        else:
            raise ValueError("spiking model not recognized")
        resp = np.max(trl.resp)
        if np.max(resp) > tp.thresh_resp:
            hi = trl.amp
        else:
            lo = trl.amp
        print(f"amp = {trl.amp} resp = {resp} thresh = {tp.thresh_resp}")

    return trl.amp


def loop_model(tp: Struct, T, nc=None):
    loop_trl = []
    for i in range(len(T)):
        row = T[i]
        trl = Struct(pw=row["pw"], amp=row["amp"], dur=row["dur"], freq=row["freq"], simdur=3)
        trl = define_trial(tp, trl)

        t = trl.t[:: tp.tSamp] if (tp.isfield("tSamp") and tp.tSamp != 1) else trl.t
        dt = t[1] - t[0]
        h = gamma(tp.ncascades, tp.tau2, t)
        csum = np.cumsum(h) * dt
        tid = np.flatnonzero(csum > 0.999)
        if tid.size:
            h = h[: tid[0] + 1]
        trl.imp_resp = h

        if tp.spikemodel == "integratefire" and nc is not None:
            result = spike_model(tp, trl, nc)
        else:
            result = spike_model(tp, trl)
        loop_trl.append(result[0] if isinstance(result, tuple) else result)
    return loop_trl


def fit_brightness(tp: Struct, T, nc=None):
    loop_trl = loop_model(tp, T, nc) if nc is not None else loop_model(tp, T)
    y_est = np.array([t.max_temporal_response for t in loop_trl])
    y = np.array([row["brightness"] for row in T])
    ind = ~np.isnan(y_est) & ~np.isnan(y)
    err = np.sum((y[ind] - y_est[ind]) ** 2)
    print(
        f"tau1 ={tp.tau1:6.4f}, tau2 ={tp.tau2:6.4f}, power ={tp.power:6.4f}, "
        f"sc_in ={tp.sc_in:6.4f}, sc_out ={tp.sc_out:6.4f}, sse = {err:6.4f}"
    )
    return err


def loop_find_threshold(tp: Struct, T, nc_or_ncells=None):
    tp.setdefault("nReps", 12)
    thresh = np.full(len(T), np.nan)
    for i, row in enumerate(T):
        trl = Struct(pw=row["pw"], amp=1, dur=row["dur"], freq=row["freq"], simdur=3)
        trl = define_trial(tp, trl)
        thresh[i] = find_threshold(trl, tp, nc_or_ncells)

    if "amp" in (T[0].keys() if T else []):
        y = np.array([row["amp"] for row in T])
        err = np.nansum((thresh - y) ** 2)
        print(f"err = {round(err, 3)}")
    else:
        err = np.nan

    if tp.isfield("experimentList"):
        for i, name in enumerate(tp.experimentList):
            print(f"{name:>10s}: {tp.sc_in[i]:g}")

    return err, thresh


# ===========================================================================
# Misc active utility (used only by external, non-Elche fitting code)
# ===========================================================================
def fit_phosphene(p: Struct, trl: Struct, v: Struct):
    """2D-Gaussian fit error (assumes circularity). Not on the Elche path."""
    phos = np.mean(trl.max_phosphene, axis=2)
    phos = phos / np.max(phos)
    pred = Gauss(p, v)
    pred = pred / np.max(pred)
    err = np.sum((phos - pred) ** 2)
    if p.sigma < 0:
        err = err + abs(p.sigma) * 1e6
    return err


def Gauss(p: Struct, v: Struct):
    return np.exp(-(((v.X - p.x) ** 2 + (v.Y - p.y) ** 2) / (2 * p.sigma**2)))
