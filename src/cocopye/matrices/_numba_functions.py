# Copyright 2023 Niklas Birth
#
# This file is part of CoCoPyE.
#
# CoCoPyE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CoCoPyE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CoCoPyE. If not, see <https://www.gnu.org/licenses/>.

"""
This module contains several numba functions that are called by functions in the parent module. In fact, these functions
mostly do the actual work while the matrix classes are just wrappers around them.
"""

from numba import njit, prange
import numpy as np
import numpy.typing as npt
from typing import Tuple, Optional


@njit(nogil=True, parallel=True)
def estimates_njit(query_mat, db_mat, k, frac_eq, progress_bar, knn_inds: npt.NDArray[np.uint64]):
    result = np.zeros((query_mat.shape[0], 3))

    for idx in prange(len(query_mat)):
        comp, cont, num = estimate_njit(db_mat, query_mat[idx], k, frac_eq, knn_inds=knn_inds[idx, :])
        result[idx] = np.array([comp, cont, num])
        progress_bar.update(1)

    return result


@njit
def estimate_njit(
        mat: npt.NDArray[np.uint8],
        vec: npt.NDArray[np.uint8],
        k: int,
        frac_eq: float = 1.0,
        knn_inds: Optional[npt.NDArray[np.uint64]] = None
) -> Tuple[float, float, int]:
    if knn_inds is None:
        knn_mat_idx = nearest_neighbors_idx_njit(mat, vec, k)[0]
        knn_mat = mat[knn_mat_idx, :]
    else:
        knn_mat = mat[knn_inds, :]

    (mode_vals, mode_nums) = mode(knn_mat)

    (mark_inds,) = np.where(
        np.logical_and(np.logical_and(mode_vals > 0, mode_vals < 255), mode_nums >= round(k * frac_eq))
    )
    n_mark = len(mark_inds)

    if n_mark == 0:
        return -1., -1., 0

    comps = np.clip(vec[mark_inds] / mode_vals[mark_inds], 0, 1)
    comp = np.mean(comps)
    cont = np.mean(vec[mark_inds] / mode_vals[mark_inds] - comps)

    return float(comp), float(cont), n_mark


@njit(parallel=True)
def nearest_neighbors_idx_njit_mat(
        db_mat: npt.NDArray[np.uint8],
        q_mat: npt.NDArray[np.uint8],
        k: int
) -> Tuple[npt.NDArray[np.int64], npt.NDArray[np.float32]]:
    knn_inds = np.zeros((q_mat.shape[0], k), dtype=np.uint64)
    knn_scores = np.zeros(q_mat.shape[0], dtype=np.float32)
    for idx in prange(q_mat.shape[0]):
        knn_inds[idx], knn_scores[idx] = nearest_neighbors_idx_njit(db_mat, q_mat[idx], k)
    return knn_inds, knn_scores


@njit(parallel=True)
def nearest_neighbors_idx_njit(
        mat: npt.NDArray[np.uint8],
        vec: npt.NDArray[np.uint8],
        k: int
) -> Tuple[npt.NDArray[np.int64], np.float32]:
    num_refs, num_count = mat.shape

    assert vec.ndim == 1, "Vector has to be 1-dimensional"
    assert vec.shape[0] == num_count, "Vector length must be equal to the number of columns of the matrix"

    eq_counts = np.zeros(num_refs)
    norm_vec = np.sqrt((mat > 0).sum(axis=1))
    for idx in prange(len(mat)):
        eq_count = np.sum(np.logical_and(np.logical_and(0 < vec, vec < 255), mat[idx] == vec))
        eq_counts[idx] = eq_count / norm_vec[idx] / np.sqrt((vec > 0).sum())

    inds = np.flip(np.argsort(eq_counts))[:k]

    return inds, eq_counts[inds].mean()


@njit
def mode(knn_mat: npt.NDArray[np.uint8]):
    mode_vals, mode_nums = np.zeros(knn_mat.shape[1], dtype=np.uint8),  np.zeros(knn_mat.shape[1], dtype=np.uint8)

    for idx, arr in enumerate(knn_mat.T):
        counts = np.bincount(arr)
        mode_vals[idx] = np.argmax(counts)
        mode_nums[idx] = np.max(counts)

    return mode_vals, mode_nums


@njit
def mean_ax0(mat):
    result = np.zeros(mat.shape[1])
    for idx, col in enumerate(mat.T):
        result[idx] = np.mean(col)

    return result


@njit
def std_ax0(mat):
    result = np.zeros(mat.shape[1])
    for idx, col in enumerate(mat.T):
        result[idx] = np.std(col)

    return result
