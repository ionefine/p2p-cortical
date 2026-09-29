"""
Python port of the MATLAB imregcode registration algorithm used by ed.m to
align simulated phosphenes with patient drawings.

Source MATLAB files ported here (all "Copyright 2023-2024 The MathWorks, Inc."):
    fftPadSize_if.m
    complexGradientImage_if.m
    fftCorrelation2D_if.m
    normalizedGradientCorrelation_if.m
    peakLocation2D_if.m
    findTranslationNGC_if.m
    logPolarResample_if.m
    findScaleRotationNGC_if.m
    resolveSimilarityRotationAmbiguityNGC_if.m

Algorithm: log-polar Fourier-Mellin registration using normalized gradient
correlation (NGC), following Tzimiropoulos, Argyriou & Stathaki, "Robust
FFT-Based Scale-Invariant Image Registration with Image Gradients," IEEE
TPAMI, vol. 32, no. 10, Oct. 2010, pp. 1899-1906.

Coordinate convention
----------------------
This module keeps MATLAB's 1-based "world == intrinsic" pixel coordinate
convention internally (pixel (row, col), 1-based, has its center at world
coordinate (x=col, y=row); an image of size (rows, cols) has default world
limits X in [0.5, cols+0.5], Y in [0.5, rows+0.5]). A forward similarity
transform is stored as a 3x3 matrix A such that, for column vectors,
    [x_out; y_out; 1] = A @ [x_in; y_in; 1]
matching MATLAB's affinetform2d/simtform2d convention. Only at the numpy
array-indexing boundary (map_coordinates) do we convert to 0-based indices.

Known limitation (documented, not hidden)
------------------------------------------
`imwarp_similarity_auto` (used internally by
`resolve_similarity_rotation_ambiguity_ngc` to warp onto an automatically
sized "tight" output canvas) reimplements MATLAB's automatic-OutputView
bounding-box logic from first principles (transform the 4 corners of the
input's world bounding box, take the axis-aligned bounding box, round the
canvas size up to the nearest whole pixel). This is the standard, unambiguous
way to do it and should match MATLAB in essentially all cases, but MATLAB's
imwarp internals for the *exact* auto-view pixel-grid snapping are not
publicly documented, so this is a principled reimplementation rather than a
byte-for-byte port. The other, fixed-canvas warp used throughout ed.py
(`imwarp_similarity_fixed`, mirroring `imwarp(img, tform, 'OutputView',
ref2d, 'FillValues', 0)`) has no such ambiguity: the output canvas there is
always the target subimage's own default reference frame.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates

__all__ = [
    "fft_pad_size",
    "complex_gradient_image",
    "fft_correlation_2d",
    "normalized_gradient_correlation",
    "peak_location_2d",
    "find_translation_ngc",
    "log_polar_resample",
    "find_scale_rotation_ngc",
    "resolve_similarity_rotation_ambiguity_ngc",
    "imwarp_similarity_fixed",
]


# ---------------------------------------------------------------------------
# fftPadSize_if
# ---------------------------------------------------------------------------
def fft_pad_size(n: int) -> int:
    """Smallest even number >= n with no prime factors greater than 7."""
    if n <= 0:
        raise ValueError("n must be a positive integer")
    n = int(n)
    if n % 2 == 1:
        n += 1
    while True:
        r = n
        for p in (2, 3, 5, 7):
            while r > 1 and r % p == 0:
                r //= p
        if r == 1:
            return n
        n += 2


# ---------------------------------------------------------------------------
# complexGradientImage_if
# ---------------------------------------------------------------------------
def complex_gradient_image(img: np.ndarray) -> np.ndarray:
    """Complex gradient image G = Gx + 1j*Gy (Tzimiropoulos 2010, eq. 12).

    Centered first-order difference in the interior, one-sided (unscaled)
    difference at the image boundary, matching MATLAB's `gradient`.
    """
    I = np.asarray(img, dtype=np.float64)
    M, N = I.shape

    Gx = np.empty_like(I)
    Gx[:, 0] = I[:, 1] - I[:, 0]
    Gx[:, -1] = I[:, -1] - I[:, -2]
    if N > 2:
        Gx[:, 1:-1] = 0.5 * (I[:, 2:] - I[:, :-2])

    Gy = np.empty_like(I)
    Gy[0, :] = I[1, :] - I[0, :]
    Gy[-1, :] = I[-1, :] - I[-2, :]
    if M > 2:
        Gy[1:-1, :] = 0.5 * (I[2:, :] - I[:-2, :])

    return Gx + 1j * Gy


# ---------------------------------------------------------------------------
# fftCorrelation2D_if
# ---------------------------------------------------------------------------
def fft_correlation_2d(A: np.ndarray, B: np.ndarray):
    """2-D correlation of A and B via FFT (B is the one being shifted).

    Returns (C, shift_x, shift_y) where shift_x[j], shift_y[i] give the
    relative shift of B with respect to A for C[i, j].
    """
    A = np.asarray(A)
    B = np.asarray(B)
    Ma, Na = A.shape
    Mb, Nb = B.shape

    shift_x = np.arange(-(Nb - 1), Na)
    shift_y = np.arange(-(Mb - 1), Ma)

    Mc = Ma + Mb - 1
    Nc = Na + Nb - 1
    Mcp = fft_pad_size(Mc)
    Ncp = fft_pad_size(Nc)

    ABConj = np.fft.fft2(A, s=(Mcp, Ncp)) * np.conj(np.fft.fft2(B, s=(Mcp, Ncp)))
    C = np.fft.ifft2(ABConj)

    ii = np.roll(np.arange(Mcp), Mb - 1)[:Mc]
    jj = np.roll(np.arange(Ncp), Nb - 1)[:Nc]
    C = C[np.ix_(ii, jj)]

    return C, shift_x, shift_y


# ---------------------------------------------------------------------------
# normalizedGradientCorrelation_if
# ---------------------------------------------------------------------------
def normalized_gradient_correlation(I1: np.ndarray, I2: np.ndarray):
    """Normalized gradient correlation of I1, I2 (Tzimiropoulos 2010, eq. 20,
    with the modified normalization noted in the MATLAB source as fixing an
    error in the paper's eq. 28)."""
    G1 = complex_gradient_image(I1)
    G2 = complex_gradient_image(I2)
    G1_bar = G1.mean()
    G2_bar = G2.mean()

    # Arguments reversed relative to the outer (I1, I2) naming, matching the
    # MATLAB source's own comment ("Reverse the order of arguments to match
    # the interface of fftCorrelation2D").
    NGC_numerator, shift_x, shift_y = fft_correlation_2d(G2 - G2_bar, G1 - G1_bar)

    NGC_denominator = np.sqrt(
        np.sum(np.abs(G1 - G1_bar) ** 2) * np.sum(np.abs(G2 - G2_bar) ** 2)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        NGC = NGC_numerator / NGC_denominator
    NGC = np.where(np.isfinite(NGC), NGC, 0)
    NGC = np.real(NGC)
    return NGC, shift_x, shift_y


# ---------------------------------------------------------------------------
# peakLocation2D_if
# ---------------------------------------------------------------------------
def peak_location_2d(F: np.ndarray):
    """Subpixel peak location of F via 2nd-order polynomial fit on the 3x3
    neighborhood around the peak. Returns 1-based (xpeak, ypeak, F_max) to
    match MATLAB's indexing convention used by the callers in this module.
    """
    F = np.asarray(F)
    M, N = F.shape

    # Column-major flatten + argmax reproduces MATLAB's max(...,'all') +
    # ind2sub tie-break (first occurrence in column-major order).
    Ffm = F.flatten(order="F")
    lin = int(np.argmax(Ffm))
    F_max = Ffm[lin]

    yi0 = lin % M
    xi0 = lin // M
    yi = yi0 + 1
    xi = xi0 + 1

    if xi == 1 or xi == N or yi == 1 or yi == M:
        return float(xi), float(yi), F_max

    u = F[yi0 - 1 : yi0 + 2, xi0 - 1 : xi0 + 2].flatten(order="F")
    x = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1], dtype=np.float64)
    y = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1], dtype=np.float64)
    X = np.column_stack([np.ones(9), x, y, x * y, x**2, y**2])

    A, *_ = np.linalg.lstsq(X, u, rcond=None)

    denom = A[3] ** 2 - 4 * A[4] * A[5]
    if denom == 0:
        return float(xi), float(yi), F_max

    x_offset = (-A[2] * A[3] + 2 * A[5] * A[1]) / denom
    y_offset = -1.0 / denom * (A[3] * A[1] - 2 * A[4] * A[2])

    x_offset = min(abs(x_offset), 0.5) * np.sign(x_offset)
    y_offset = min(abs(y_offset), 0.5) * np.sign(y_offset)

    return xi + x_offset, yi + y_offset, F_max


# ---------------------------------------------------------------------------
# findTranslationNGC_if
# ---------------------------------------------------------------------------
def find_translation_ngc(moving: np.ndarray, fixed: np.ndarray):
    """Find the translation that shifts `moving` to align it with `fixed`
    using normalized gradient correlation. Returns ((dx, dy), peak)."""
    NGC, shift_x, shift_y = normalized_gradient_correlation(moving, fixed)
    xpeak, ypeak, peak = peak_location_2d(NGC)
    peak = float(peak)

    dx = float(np.interp(xpeak, np.arange(1, len(shift_x) + 1), shift_x))
    dy = float(np.interp(ypeak, np.arange(1, len(shift_y) + 1), shift_y))

    return (dx, dy), peak


# ---------------------------------------------------------------------------
# logPolarResample_if
# ---------------------------------------------------------------------------
def log_polar_resample(F: np.ndarray):
    """Log-polar resampling of the (real, non-negative) FFT-magnitude image F.

    Returns (L, b, thetad): L is the log-polar-resampled image, b is the
    exponential base relating a horizontal shift to a scale factor
    (s = b**shift_x), and thetad gives the angle (degrees) for each row.
    """
    F = np.asarray(F, dtype=np.float64)
    N = F.shape[0]  # F assumed N x N, N even

    K = max(N // 2, 256)
    rho = np.arange(K)
    b = (N / 2) ** (1.0 / (K - 1))

    thetad = np.arange(1800) / 10.0  # 0, 0.1, ..., 179.9

    Fs = np.fft.fftshift(F)
    x_c = N / 2 + 1
    y_c = N / 2 + 1

    r = b**rho
    theta_rad = np.deg2rad(thetad)
    xq = np.outer(np.cos(theta_rad), r) + x_c
    yq = np.outer(np.sin(theta_rad), r) + y_c

    xq = np.clip(xq, 1, N)
    yq = np.clip(yq, 1, N)

    coords = np.stack([yq - 1.0, xq - 1.0], axis=0)  # 0-based (row, col)
    L = map_coordinates(Fs, coords, order=1, mode="constant", cval=0.0)

    return L, b, thetad


# ---------------------------------------------------------------------------
# findScaleRotationNGC_if
# ---------------------------------------------------------------------------
def find_scale_rotation_ngc(I1: np.ndarray, I2: np.ndarray):
    """Find the scale factor s and rotation angle r (degrees) to align image
    I1 with image I2, using log-polar normalized gradient correlation."""
    I1 = np.asarray(I1, dtype=np.float64)
    I2 = np.asarray(I2, dtype=np.float64)

    Y1, X1 = I1.shape
    Y2, X2 = I2.shape
    Np = max(X1, Y1, X2, Y2)
    N = fft_pad_size(Np)

    F1 = np.abs(np.fft.fft2(complex_gradient_image(I1), s=(N, N)))
    F2 = np.abs(np.fft.fft2(complex_gradient_image(I2), s=(N, N)))

    L1, base, thetad = log_polar_resample(F1)
    L2, _, _ = log_polar_resample(F2)

    (tx, ty), _peak = find_translation_ngc(L1, L2)

    s = base ** (-tx)
    P = len(thetad)
    if ty >= 0:
        r = np.interp(ty, np.arange(P), thetad)
    else:
        r = -np.interp(-ty, np.arange(P), thetad)

    return float(s), float(r)


# ---------------------------------------------------------------------------
# imwarp equivalents (MATLAB affinetform2d + imwarp, similarity-only)
# ---------------------------------------------------------------------------
def _similarity_matrix(S: float, theta_deg: float, tx: float = 0.0, ty: float = 0.0):
    c = np.cos(np.deg2rad(theta_deg))
    s = np.sin(np.deg2rad(theta_deg))
    return np.array(
        [
            [S * c, -S * s, tx],
            [S * s, S * c, ty],
            [0.0, 0.0, 1.0],
        ]
    )


def _transform_corners(A: np.ndarray, cols: int, rows: int):
    """Forward-transform the 4 corners of an image's default world bounding
    box [0.5, cols+0.5] x [0.5, rows+0.5]."""
    xw = np.array([0.5, cols + 0.5, cols + 0.5, 0.5])
    yw = np.array([0.5, 0.5, rows + 0.5, rows + 0.5])
    pts = np.vstack([xw, yw, np.ones(4)])
    out = A @ pts
    return out[0], out[1]


def _warp_to_canvas(
    moving: np.ndarray,
    A: np.ndarray,
    out_rows: int,
    out_cols: int,
    x_world_limits,
    y_world_limits,
    fill_value: float = 0.0,
):
    """Resample `moving` (default world==intrinsic reference) onto an output
    canvas of size (out_rows, out_cols) whose world limits are given, using
    the inverse of forward transform A (bilinear interpolation)."""
    Ainv = np.linalg.inv(A)

    col_idx = np.arange(1, out_cols + 1)
    row_idx = np.arange(1, out_rows + 1)
    x_world = x_world_limits[0] + (col_idx - 0.5)
    y_world = y_world_limits[0] + (row_idx - 0.5)
    Xw, Yw = np.meshgrid(x_world, y_world)

    pts = np.vstack([Xw.ravel(), Yw.ravel(), np.ones(Xw.size)])
    src = Ainv @ pts
    x_in = src[0].reshape(Xw.shape)
    y_in = src[1].reshape(Xw.shape)

    col_in = x_in - 1.0  # moving image: default world == intrinsic
    row_in = y_in - 1.0
    coords = np.stack([row_in, col_in], axis=0)

    return map_coordinates(
        np.asarray(moving, dtype=np.float64),
        coords,
        order=1,
        mode="constant",
        cval=fill_value,
    )


def imwarp_similarity_auto(moving: np.ndarray, A: np.ndarray):
    """Equivalent of MATLAB `imwarp(moving, affinetform2d(A), 'SmoothEdges', true)`
    with no OutputView specified: the output canvas is automatically sized to
    tightly bound the transformed image (see module docstring for the
    documented limitation on exact pixel-grid-snapping fidelity).

    Returns (warped_image, x_world_limits, y_world_limits).
    """
    rows, cols = moving.shape
    xo, yo = _transform_corners(A, cols, rows)
    x_lo, x_hi = float(xo.min()), float(xo.max())
    y_lo, y_hi = float(yo.min()), float(yo.max())

    out_cols = max(int(np.ceil(x_hi - x_lo)), 1)
    out_rows = max(int(np.ceil(y_hi - y_lo)), 1)

    x_world_limits = (x_lo, x_lo + out_cols)
    y_world_limits = (y_lo, y_lo + out_rows)

    warped = _warp_to_canvas(
        moving, A, out_rows, out_cols, x_world_limits, y_world_limits
    )
    return warped, x_world_limits, y_world_limits


def imwarp_similarity_fixed(
    moving: np.ndarray, A: np.ndarray, out_shape: tuple[int, int], fill_value: float = 0.0
):
    """Equivalent of MATLAB `imwarp(moving, tform, 'OutputView', imref2d(out_shape),
    'FillValues', fill_value)`: warp onto a fixed-size output canvas using its
    default (world == intrinsic) reference frame. `out_shape` is (rows, cols),
    matching MATLAB's `size(target)` / `imref2d(targetSize)` convention.
    """
    out_rows, out_cols = out_shape
    x_world_limits = (0.5, out_cols + 0.5)
    y_world_limits = (0.5, out_rows + 0.5)
    return _warp_to_canvas(
        moving, A, out_rows, out_cols, x_world_limits, y_world_limits, fill_value
    )


# ---------------------------------------------------------------------------
# resolveSimilarityRotationAmbiguityNGC_if
# ---------------------------------------------------------------------------
def resolve_similarity_rotation_ambiguity_ngc(
    moving: np.ndarray, fixed: np.ndarray, S: float, thetad: float
):
    """Find the best translation and final rotation angle to align `moving`
    to `fixed` given a scale and candidate rotation angle that might be off
    by 180 degrees. Returns (A, peak) where A is a 3x3 forward similarity
    transform matrix (see module docstring for convention)."""
    moving = np.asarray(moving, dtype=np.float64)
    fixed = np.asarray(fixed, dtype=np.float64)

    flip_inputs = S < 1
    if flip_inputs:
        S = 1.0 / S
        thetad = -thetad
        moving, fixed = fixed, moving  # matches the MATLAB swap() trick

    theta1d = thetad
    theta2d = thetad + 180.0

    A1 = _similarity_matrix(S, theta1d)
    A2 = _similarity_matrix(S, theta2d)

    scaledRotatedMoving1, xwl1, ywl1 = imwarp_similarity_auto(moving, A1)
    scaledRotatedMoving2 = np.rot90(scaledRotatedMoving1, 2)
    xwl2 = tuple(sorted((-xwl1[0], -xwl1[1])))
    ywl2 = tuple(sorted((-ywl1[0], -ywl1[1])))

    (vec1x, vec1y), peak1 = find_translation_ngc(scaledRotatedMoving1, fixed)
    (vec2x, vec2y), peak2 = find_translation_ngc(scaledRotatedMoving2, fixed)

    if peak1 >= peak2:
        vec = (vec1x, vec1y)
        A = A1.copy()
        xwl, ywl = xwl1, ywl1
        peak = peak1
    else:
        vec = (vec2x, vec2y)
        A = A2.copy()
        xwl, ywl = xwl2, ywl2
        peak = peak2

    # XIntrinsicLimits(1)/YIntrinsicLimits(1) are always 0.5 for any imref2d canvas.
    finalXOffset = vec[0] + (0.5 - xwl[0])
    finalYOffset = vec[1] + (0.5 - ywl[0])

    A[0, 2] = finalXOffset
    A[1, 2] = finalYOffset

    if flip_inputs:
        A = np.linalg.inv(A)

    return A, peak
