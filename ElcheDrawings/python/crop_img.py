"""
Python port of crop_img.m (Guilherme Coco Beltramini, 2013-May-29).

Crops an image by removing homogeneous-intensity borders (background),
optionally leaving a margin of `border` pixels.
"""

from __future__ import annotations

import numpy as np

__all__ = ["crop_img"]


def crop_img(img: np.ndarray, border: int = 0):
    """Crop image by removing edges with homogeneous intensity.

    Parameters
    ----------
    img : ndarray
        MxN or MxNxC array (C = number of color layers).
    border : int
        Maximum number of pixels to leave at the borders (default 0).

    Returns
    -------
    img2 : ndarray
        Cropped image, same dtype and number of dims as `img`.
    edge_row1, edge_col1, edge_row2, edge_col2 : int
        1-based (MATLAB-convention) row/col bounds of the crop, for parity
        with the MATLAB function's extra return values.
    """
    img = np.asarray(img)
    squeeze_back = img.ndim == 2
    if squeeze_back:
        img = img[:, :, np.newaxis]

    MM, NN, CC = img.shape
    edge_col = np.zeros((2, CC), dtype=np.int64)
    edge_row = np.zeros((2, CC), dtype=np.int64)

    for cc in range(CC):
        layer = img[:, :, cc]

        # --- Top-left corner ---
        img_bg = layer == layer[0, 0]

        cols = img_bg.sum(axis=0)  # length NN
        if cols[0] == MM:
            d = np.diff(cols.astype(np.int64))
            nz = np.flatnonzero(d)
            if nz.size:
                tmp = nz[0] + 1  # 1-based MATLAB "find(diff(cols),1,'first')"
                edge_col[0, cc] = tmp + 1 - border
            else:
                edge_col[0, cc] = 1
        else:
            edge_col[0, cc] = 1

        rows = img_bg.sum(axis=1)  # length MM
        if rows[0] == NN:
            d = np.diff(rows.astype(np.int64))
            nz = np.flatnonzero(d)
            if nz.size:
                tmp = nz[0] + 1
                edge_row[0, cc] = tmp + 1 - border
            else:
                edge_row[0, cc] = 1
        else:
            edge_row[0, cc] = 1

        # --- Bottom-right corner ---
        img_bg = layer == layer[MM - 1, NN - 1]

        cols = img_bg.sum(axis=0)
        if cols[-1] == MM:
            d = np.diff(cols.astype(np.int64))
            nz = np.flatnonzero(d)
            if nz.size:
                tmp = nz[-1] + 1  # "find(diff(cols),1,'last')", 1-based
                edge_col[1, cc] = tmp + border
            else:
                edge_col[1, cc] = NN
        else:
            edge_col[1, cc] = NN

        rows = img_bg.sum(axis=1)
        if rows[-1] == NN:
            d = np.diff(rows.astype(np.int64))
            nz = np.flatnonzero(d)
            if nz.size:
                tmp = nz[-1] + 1
                edge_row[1, cc] = tmp + border
            else:
                edge_row[1, cc] = MM
        else:
            edge_row[1, cc] = MM

        # --- Identify homogeneous color layers (ignore them) ---
        if (
            edge_col[0, cc] == 1
            and edge_col[1, cc] == NN
            and edge_row[0, cc] == 1
            and edge_row[1, cc] == MM
            and not np.any(np.diff(layer.astype(np.int64), axis=0))
            and not np.any(np.diff(layer.astype(np.int64), axis=1))
        ):
            edge_col[:, cc] = [NN, 1]
            edge_row[:, cc] = [MM, 1]

    # --- Combine across color layers ---
    edge_col_lo = edge_col[0, :].min()
    edge_col_hi = edge_col[1, :].max()
    edge_col = np.array([edge_col_lo, edge_col_hi])
    edge_col[0] = max(edge_col[0], 1)
    edge_col[1] = min(edge_col[1], NN)

    edge_row_lo = edge_row[0, :].min()
    edge_row_hi = edge_row[1, :].max()
    edge_row = np.array([edge_row_lo, edge_row_hi])
    edge_row[0] = max(edge_row[0], 1)
    edge_row[1] = min(edge_row[1], MM)

    # --- Crop (convert 1-based inclusive MATLAB bounds to 0-based numpy slice) ---
    r1, r2 = int(edge_row[0]), int(edge_row[1])
    c1, c2 = int(edge_col[0]), int(edge_col[1])
    img2 = img[r1 - 1 : r2, c1 - 1 : c2, :]

    if squeeze_back:
        img2 = img2[:, :, 0]

    return img2, r1, c1, r2, c2
