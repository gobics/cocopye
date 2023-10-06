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
        frac_eq: float = 0.9,
        knn_inds: Optional[npt.NDArray[np.uint64]] = None
) -> Tuple[float, float, int]:
    """
    Calculate an estimate for completeness and contamination for a vector based on common markers in the k nearest
    neighbors.
    :param mat: Database matrix that will be used for the nearest neighbors.
    :param vec: Input vector
    :param k: k (like in k nearest neighbors). This parameter is ignored if knn_inds is provided.
    :param frac_eq: Fraction of similar counts within the nearest neighbors required to consider a Pfam/kmer as a
    marker
    :param knn_inds: If provided, this is used as the k nearest neighbors (1D-array containing indices of mat).
    :return: A 3-tuple: First element is the completeness estimate, the second is the contamiation estimate
    (both between 0 and 1) and the third one is the number of markers that were used. This last value is mainly
    intended for evaluation purposes.
    """
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
def nearest_neighbors_idx_njit(mat: npt.NDArray[np.uint8], vec: npt.NDArray[np.uint8], k: int) -> Tuple[npt.NDArray[np.int64], np.float32]:
    """
    Returns the row indices of the k nearest neighbors of a vector in the databse matrix. This is mainly used by the
    `nearest_neighbors`  function, but may also be useful in other situation where one needs only the indices and
    not the actual neighbors.

    Parameters are the same as of `nearest_neighbors`.

    :return: A numpy array containing the indices of the nearest neighbors.
    """
    num_refs, num_count = mat.shape

    assert vec.ndim == 1, "Vector has to be 1-dimensional"
    assert vec.shape[0] == num_count, "Vector length must be equal to the number of columns of the matrix"

    eq_count = np.zeros(num_refs)
    norm_vec = np.sqrt((mat > 0).sum(axis=1))
    for idx in prange(len(mat)):
        eq_count[idx] = np.sum(np.logical_and(np.logical_and(0 < vec, vec < 255), mat[idx] == vec)) / norm_vec[idx] / np.sqrt((vec > 0).sum())

    inds = np.flip(np.argsort(eq_count))[:k]

    return inds, eq_count[inds].mean()


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
