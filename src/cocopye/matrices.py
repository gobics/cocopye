"""
This module contains several matrix classes that are relevant for our work. They are based on numpy arrays.

**Matrix** is a generic class and the ancestor of the other matrix classes. You probably don't need to use it directly.

**DatabaseMatrix** and **QueryMatrix** are quite similar in the sense that they both represent Pfam or Kmer count
matrices. However, for better differentiation I decided to split them up into two classes.

A **DatabaseMatrix**, as the name implies, contains the count values of our database sequences, while a
**QueryMatrix** contains the counts of the input bins (the one we want to determine completeness and contamination for).
"""
from __future__ import annotations

import pickle
from typing import TypeVar, Generic, cast, Tuple, Optional
import numpy as np
import numpy.typing as npt
from numba import njit, prange
from numba_progress import ProgressBar

from .histogram import Histogram

T = TypeVar("T")


class Matrix(Generic[T]):
    """
    A generic class for a 2-dimensional matrix. It is intended as an ancestor for the other matrix classes in this
    module.
    """
    _mat: T

    def __init__(self, mat: T):
        cast(npt.NDArray[np.generic], mat)
        assert mat.ndim == 2, "Matrix has to be 2-dimensional"  # type: ignore

        self._mat = mat

    def mat(self) -> T:
        """
        This returns the underlying numpy array in case direct access to it is required.
        """
        return self._mat

    def save_to_file(self, filename: str) -> None:
        """
        Save the matrix to a file. If the filename has a .npy extension it is saved in binary format, otherwise as csv.
        This is just a convenience function. If it doesn't fit to your needs, just use np.save or np.savetxt.
        """
        file_format = filename.split(".")[-1]

        if file_format == "npy":
            np.save(filename, self._mat)  # type: ignore
        else:
            np.savetxt(filename, self._mat, delimiter=",")  # type: ignore

    def __str__(self) -> str:
        return str(self._mat)


class DatabaseMatrix(Matrix[npt.NDArray[np.uint8]]):
    """
    Description
    """
    def __init__(self, mat: npt.NDArray[np.uint8]):
        """
        :param mat: 2-dimensional numpy matrix. Rows are expected to represent sequences/bins while columns are Pfams or
        kmers.
        """
        super().__init__(mat)

    def nearest_neighbors(self, vec: npt.NDArray[np.uint8], k: int) -> DatabaseMatrix:
        """
        Returns the k nearest neighbors of a vector in the database matrix.
        :param vec: The count vector you want to find neighbors for
        :param k: Well, it's the k of k nearest neighbors
        :return: Database matrix where each row is one of the nearest neighbors.
        """
        return DatabaseMatrix(self._mat[self.nearest_neighbors_idx(vec, k), :])

    def nearest_neighbors_idx(self, vec: npt.NDArray[np.uint8], k: int) -> npt.NDArray[np.int64]:
        """
        Returns the row indices of the k nearest neighbors of a vector in the databse matrix. This is mainly used by the
        `nearest_neighbors`  function, but may also be useful in other situation where one needs only the indices and
        not the actual neighbors.

        Parameters are the same as of `nearest_neighbors`.

        :return: A numpy array containing the indices of the nearest neighbors.
        """
        return nearest_neighbors_idx_njit(self._mat, vec, k)

    def universal_markers(self, threshold: float = 0.95) -> npt.NDArray[np.uint32]:
        return np.where(np.count_nonzero(self._mat == 1, axis=0) / self._mat.shape[0] >= threshold)[0]

    def estimate(
            self,
            vec: npt.NDArray[np.uint8],
            k: int,
            frac_eq: float = 0.9,
            knn: Optional[DatabaseMatrix] = None,
    ) -> Tuple[float, float, int]:
        """
        Calculate an estimate for completeness and contamination for a vector based on common markers in the k nearest
        neighbors.
        :param vec: Input vector
        :param k: k (like in k nearest neighbors)
        :param frac_eq: Fraction of similar counts within the nearest neighbors required to consider a Pfam/kmer as a
        marker
        :param knn: If provided, this is used as the k nearest neighbors
        :return: A 3-tuple: First element is the completeness estimate, the second is the contamiation estimate
        (both between 0 and 1) and the third one is the number of markers that were used. This last value is mainly
        intended for evaluation purposes.
        """
        if knn is not None:
            knn = knn.mat()
        return estimate_njit(self.mat(), vec, k, frac_eq, knn)


class QueryMatrix(Matrix[npt.NDArray[np.uint8]]):
    """
    Description
    """
    def __init__(self, mat: npt.NDArray[np.uint8]):
        """
        :param mat: 2-dimensional numpy matrix. Rows are expected to represent sequences/bins while columns are Pfams or
        kmers.
        """
        super().__init__(mat)

    def knn_inds(self, db: DatabaseMatrix, k: int):
        return nearest_neighbors_idx_njit_mat(db.mat(), self._mat, k)

    def preestimates(self, markers: npt.NDArray[np.uint32]) -> npt.NDArray[np.flot32]:
        submatrix = self._mat[:, markers]

        completeness = np.sum(np.clip(submatrix, 0, 1), axis=1) / markers.shape[0]
        contamination = np.sum(submatrix - np.clip(submatrix, 0, 1), axis=1) / markers.shape[0]

        return np.array([completeness, contamination], dtype=np.float32).T

    def estimates(self, db: DatabaseMatrix, k: int, knn_inds: npt.NDArray[np.uint64], frac_eq: float = 0.9) -> npt.NDArray[np.float32]:
        """
        Calculate a completeness and contamination estimate for all rows in the QueryMatrix based on common markers in
        the k nearest neighbors.
        :param db: Database used to obtain the neighbors
        :param k: k for k nearest neighbors
        :param frac_eq: Fraction of similar counts within the nearest neighbors required to consider a Pfam/kmer as a
        marker
        :return: A 2-dimensional numpy array. For each row in the QueryMatrix there is a row with two floats, where the
        first element ist the completeness and the second one the contamination estimate.
        """
        with ProgressBar(total=self.mat().shape[0], ncols=100, dynamic_ncols=False, desc="- Calculating estimates") as progress_bar:
            result = estimates_njit(self.mat(), db.mat(), k, frac_eq, progress_bar, knn_inds)
        return result

    def into_feature_mat(self, db: DatabaseMatrix, estimates: npt.NDArray[np.float32], knn_inds: npt.NDArray[np.uint64], resolution: int = 10) -> FeatureMatrix:
        hist = Histogram(resolution)

        feature_vecs = []
        for idx, row in enumerate(self._mat):
            knn_mat = db.mat()[knn_inds[idx]]

            knn_hists = []
            for neighbor in knn_mat:
                knn_hists.append(hist.calc_bins_for_two_counts(row, neighbor))

            mean_hist = np.mean(np.array(knn_hists), axis=0)
            mean_hist = mean_hist / np.sum(mean_hist)

            feature_vecs.append(np.concatenate([mean_hist, estimates[idx, :2]]))

        return FeatureMatrix(np.array(feature_vecs))


class FeatureMatrix(Matrix[npt.NDArray[np.double]]):
    def __init__(self, mat: npt.NDArray[np.double]):
        super().__init__(mat)

    def ml_estimates(self, model_file: str):
        mlp = pickle.load(open(model_file, "rb"))
        return mlp.predict(self._mat)


def load_u8mat_from_file(filename: str) -> npt.NDArray[np.uint8]:
    """
    This is just a convenience function to load a numpy matrix from a file. If it doesn't suit your requirements, just
    use numpy.load or numpy.loadtxt directly.
    :param filename: Filename of the matrix file. If the extension is .npy it is assumed that the content is in binary
    format, otherwise it should be csv.
    :return: The loaded matrix as a numpy array
    """
    file_format = filename.split(".")[-1]

    mat: npt.NDArray[np.uint8]
    if file_format == "npy":
        mat = np.load(filename)
    else:
        mat = np.loadtxt(filename, delimiter=",", dtype=np.uint8)

    return mat


# === NJIT FUNCTIONS ===================================================================================================

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
    :param vec: Input vector
    :param k: k (like in k nearest neighbors)
    :param frac_eq: Fraction of similar counts within the nearest neighbors required to consider a Pfam/kmer as a
    marker
    :param knn: If provided, this is used as the k nearest neighbors
    :return: A 3-tuple: First element is the completeness estimate, the second is the contamiation estimate
    (both between 0 and 1) and the third one is the number of markers that were used. This last value is mainly
    intended for evaluation purposes.
    """
    if knn_inds is None:
        knn_mat_idx = nearest_neighbors_idx_njit(mat, vec, k)
        knn_mat = mat[knn_mat_idx, :]
    else:
        knn_mat = mat[knn_inds, :]

    (mode_vals, mode_nums) = mode(knn_mat)

    (mark_inds,) = np.where(np.logical_and(np.logical_and(mode_vals > 0, mode_vals < 255), mode_nums >= round(k * frac_eq)))
    n_mark = len(mark_inds)

    if n_mark == 0:
        return -1., -1., 0

    comps = np.clip(vec[mark_inds] / mode_vals[mark_inds], 0, 1)
    comp = np.mean(comps)
    cont = np.mean(vec[mark_inds] / mode_vals[mark_inds] - comps)

    return float(comp), float(cont), n_mark


@njit(parallel=True)
def nearest_neighbors_idx_njit_mat(db_mat: npt.NDArray[np.int8], q_mat: npt.NDArray[np.uint8], k: int) -> npt.NDArray[np.int64]:
    knn_inds = np.zeros((q_mat.shape[0], k), dtype=np.uint64)
    for idx in prange(q_mat.shape[0]):
        knn_inds[idx] = nearest_neighbors_idx_njit(db_mat, q_mat[idx], k)
    return knn_inds


@njit(parallel=True)
def nearest_neighbors_idx_njit(mat: npt.NDArray[np.int8], vec: npt.NDArray[np.uint8], k: int) -> npt.NDArray[np.int64]:
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
    for idx in prange(len(mat)):
        eq_count[idx] = np.sum(np.logical_and((0 < vec) < 255, mat[idx] == vec))

    return np.flip(np.argsort(eq_count))[:k]


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


# ======================================================================================================================
#
# The following contains some functions that are no longer (or not yet) used. We probably can delete them in the future.
#
# Functions of DatabaseMatrix
#
# def universal_markers_idx(self, threshold: float) -> npt.NDArray[np.uint64]:
#     return np.where(np.count_nonzero(self._mat == 1, axis=0) / self._mat.shape[0] >= threshold)[0]

# def preestimate(self, markers: npt.NDArray[np.uint64], vec: npt.NDArray[np.uint8]) -> (float, float):
#     # TODO: Wie soll das genau berechnet werden? Completeness nur bei marker == 1 oder marker >= 1?
#     # Contamination bei marker > 1? Oder gar keine Vorabschätzung für Contamination?
#     pass
#
# Functions of QueryMatrix
#
# def preestimates(self, db: DatabaseMatrix, threshold: float) -> npt.NDArray[np.float]:
#     pass
#
# def transform_features(self, vec: npt.NDArray[np.uint8], n: int = 100, minmax: float = 3.) -> FeatureMatrix:
#     eps = 1e-8
#     n_rows = self._mat.shape[0]
#     transformed_mat = np.zeros((n_rows, n))
#
#    for i in range(n_rows):
#         # Normal values (where at least one part is > 0 and we have no 255)
#         inds = np.where(np.logical_and(np.logical_or(vec > 0, self._mat[i] > 0), vec < 255 and self._mat[i] < 255))[0]
#         normal_values = np.log2((vec[inds] + eps) / (self._mat[i, inds] + eps))
#
#         # Special values where at least one part contains max value (255)
#         inds = np.where(np.logical_or(vec == 255, self._mat[i] == 255))[0]
#         max_values = np.log2((vec[inds] + eps) / (self._mat[i, inds] + eps))
#
#         # Count all the values and store them in our matrix
#         transformed_mat[i, 0] += np.sum(np.concatenate((normal_values, max_values)) > minmax)
#         transformed_mat[i, 1:] = np.histogram(normal_values, n - 2, range=(-minmax, minmax))[0]
#         transformed_mat[i, -1] += np.sum(np.concatenate((normal_values, max_values)) < -minmax)
#
#         # Do we need this? TODO
#         # transformed_mat[i] /= np.sum(transformed_mat[i])
#
#     return FeatureMatrix(transformed_mat)
#
# Class FeatureMatrix
#
# class FeatureMatrix(Matrix[npt.NDArray[np.double]]):
#     def __init__(self, mat: npt.NDArray[np.double]):
#         super().__init__(mat)
#
#     def neural_net_estimates(self):
#         # Ist diese Funktion so sinnvoll?
#         pass
