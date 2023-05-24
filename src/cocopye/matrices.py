"""
This module contains several matrix classes that are relevant for our work. They are based on numpy arrays.

**Matrix** is a generic class and the ancestor of the other matrix classes. You probably don't need to use it directly.

**DatabaseMatrix** and **QueryMatrix** are quite similar in the sense that they both represent Pfam or Kmer count
matrices. However, for better differentiation I decided to split them up into two classes.

A **DatabaseMatrix**, as the name implies, contains the count values of our database sequences, while a
**QueryMatrix** contains the counts of the input bins (the one we want to determine completeness and contamination for).
"""
from __future__ import annotations
from typing import TypeVar, Generic, cast, Tuple, Optional
import numpy as np
import numpy.typing as npt
import scipy.stats as st


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
    TODO: Description
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
        num_refs, num_count = self._mat.shape

        assert vec.ndim == 1, "Vector has to be 1-dimensional"
        assert vec.shape[0] == num_count, "Vector length must be equal to the number of columns of the matrix"

        eq_count = np.zeros(num_refs)
        for idx, row in enumerate(self._mat):
            eq_count[idx] = np.sum(np.logical_and((0 < vec) < 255, row == vec))

        return np.flip(np.argsort(eq_count))[:k]

    def estimate(
            self,
            vec: npt.NDArray[np.uint8],
            k: int,
            frac_eq: float = 0.9,
            knn: Optional[DatabaseMatrix] = None
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
        knn_mat = self.nearest_neighbors(vec, k).mat() if knn is None else knn.mat()
        (mode_vals, mode_nums) = st.mode(knn_mat, axis=0, keepdims=False)

        (mark_inds,) = np.where(np.logical_and(mode_vals > 0, mode_nums >= round(k * frac_eq)))
        n_mark = len(mark_inds)

        comps = np.clip(vec[mark_inds] / mode_vals[mark_inds], 0, 1)
        comp = np.mean(comps)
        cont = np.mean(vec[mark_inds] / mode_vals[mark_inds] - comps)

        return float(comp), float(cont), n_mark


class QueryMatrix(Matrix[npt.NDArray[np.uint8]]):
    """
    TODO: Description
    """
    def __init__(self, mat: npt.NDArray[np.uint8]):
        """
        :param mat: 2-dimensional numpy matrix. Rows are expected to represent sequences/bins while columns are Pfams or
        kmers.
        """
        super().__init__(mat)

    def estimates(self, db: DatabaseMatrix, k: int, frac_eq: float = 0.9) -> npt.NDArray[np.float32]:
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
        def func(vec: npt.NDArray[np.uint8], k_inner: int, frac_eq_inner: float) -> npt.NDArray[np.float32]:
            comp, cont, num = db.estimate(vec, k_inner, frac_eq_inner)
            return np.array([comp, cont])

        return np.apply_along_axis(func, 1, self.mat(), k, frac_eq)


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