"""Fixed-radius neighbour search.

Uses scipy's cKDTree when available and falls back to a numpy grid hash, so
the plugin works in KiCad's bundled Python where only numpy may exist.
"""
from __future__ import annotations

import itertools

import numpy as np

try:  # pragma: no cover - exercised implicitly
    from scipy.spatial import cKDTree as _KD
except Exception:  # pragma: no cover
    _KD = None


def query_pairs(points: np.ndarray, r: float, force_numpy: bool = False) -> np.ndarray:
    """All index pairs (i < j) with |p_i - p_j| <= r, as an (M, 2) int array."""
    pts = np.asarray(points, float)
    if len(pts) < 2:
        return np.zeros((0, 2), int)
    if _KD is not None and not force_numpy:
        return _KD(pts).query_pairs(r, output_type="ndarray")
    d = pts.shape[1]
    cell = np.floor((pts - pts.min(axis=0)) / r).astype(np.int64)
    dims = cell.max(axis=0) + 3
    mult = np.cumprod(np.concatenate([[1], dims[:-1]]))
    key = (cell + 1) @ mult
    order = np.argsort(key, kind="stable")
    ks = key[order]
    uniq, start = np.unique(ks, return_index=True)
    end = np.append(start[1:], len(ks))
    lookup = dict(zip(uniq.tolist(), zip(start.tolist(), end.tolist())))
    offs = [np.array(o) for o in itertools.product((-1, 0, 1), repeat=d)]
    offs = [o for o in offs if tuple(o) >= tuple([0] * d)]      # half space incl. zero
    out = []
    r2 = r * r
    for kc, (a, b) in lookup.items():
        ia = order[a:b]
        for o in offs:
            kn = kc + int(o @ mult)
            if kn not in lookup:
                continue
            c, e = lookup[kn]
            ib = order[c:e]
            diff = pts[ia][:, None, :] - pts[ib][None, :, :]
            dd = np.einsum("ijk,ijk->ij", diff, diff)
            ii, jj = np.nonzero(dd <= r2)
            gi, gj = ia[ii], ib[jj]
            if not o.any():
                m = gi < gj
                gi, gj = gi[m], gj[m]
            out.append(np.stack([np.minimum(gi, gj), np.maximum(gi, gj)], axis=1))
    if not out:
        return np.zeros((0, 2), int)
    return np.concatenate(out)
