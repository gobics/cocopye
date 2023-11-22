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
This module contains several matrix classes that are relevant for our work. They are based on numpy arrays.

**Matrix** is a generic class and the ancestor of the other matrix classes. You probably don't need to use it directly.

**DatabaseMatrix** and **QueryMatrix** are quite similar in the sense that they both represent Pfam count
matrices. However, for better differentiation I decided to split them up into two classes.

A **DatabaseMatrix**, as the name implies, contains the count values of our database sequences, while a
**QueryMatrix** contains the counts of the input bins (the one we want to determine completeness and contamination for).

A DatabaseMatrix can be transformed into a **FeatureMatrix** to get completeness and contamination estimates from the
neural network.
"""

from __future__ import annotations

import pickle
from datetime import datetime
from typing import TypeVar, Generic, cast, Tuple, Optional, List
import numpy as np
import numpy.typing as npt
from numba_progress import ProgressBar
import pandas as pd

from ..histogram import Histogram
from ._numba_functions import nearest_neighbors_idx_njit, nearest_neighbors_idx_njit_mat, estimates_njit, estimate_njit

T = TypeVar("T")


class Matrix(Generic[T]):
    """
    A generic class for a 2-dimensional matrix. It is intended as an ancestor for the other matrix classes in this
    module.

    Please note that even though we are using a completely generic parameter, the class must be used with numpy arrays.
    """
    _mat: T

    def __init__(self, mat: T):
        """
        :param mat: Some kind of two-dimensional numpy matrix
        """
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
        Save the matrix to a file. If the filename has a .npy extension it is saved in binary format. If the extension
        is .npz it uses some form of compression. Otherwise the matrix is stored as csv.

        This is just a convenience function. If it doesn't suit your needs, just use np.save or np.savetxt.
        """
        file_format = filename.split(".")[-1]

        if file_format == "npy":
            np.save(filename, self._mat)  # type: ignore
        elif file_format == "npz":
            np.savez_compressed(filename, self._mat)
        else:
            np.savetxt(filename, self._mat, delimiter=",")  # type: ignore

    def __str__(self) -> str:
        return str(self._mat)


class DatabaseMatrix(Matrix[npt.NDArray[np.uint8]]):
    """
    A two-dimensional matrix containing Pfam counts. It can optionally contain metadata about the underlying sequences.
    """
    _metadata: Optional[pd.DataFrame] = None

    def __init__(self, mat: npt.NDArray[np.uint8], metadata: Optional[pd.DataFrame] = None):
        """
        :param mat: 2-dimensional numpy matrix. Each row contains the counts of a reference sequences while each column
        represents a Pfam.
        :param metadata: A Pandas dataframe containing metadata for each reference sequence. It is expected that the
        same row indices of the matrix and metadata represent the same sequence.
        """
        super().__init__(mat)
        self._metadata = metadata

    def metadata(self) -> Optional[pd.DataFrame]:
        """
        :return: The dataframe with metadata (in case you need direct access to it). This will return None if the
        DatabaseMatrix contains no metadata.
        """
        return self._metadata

    def nearest_neighbors(self, vec: npt.NDArray[np.uint8], k: int) -> DatabaseMatrix:
        """
        Returns the k nearest neighbors of a vector in the database matrix.
        :param vec: The count vector you want to find neighbors for
        :param k: Number of nearest neighbors to return
        :return: Database matrix where each row is one of the nearest neighbors.
        """
        return DatabaseMatrix(self._mat[self.nearest_neighbors_idx(vec, k), :])

    def nearest_neighbors_idx(self, vec: npt.NDArray[np.uint8], k: int) -> npt.NDArray[np.int64]:
        """
        Returns the row indices of the k nearest neighbors of a vector in the databse matrix. This is mainly used by the
        `nearest_neighbors`  function, but may also be useful in other situations where one needs only the indices and
        not the actual neighbors.

        :param vec: The count vector you want to find neighbors for
        :param k: Number of nearest neighbors to return

        :return: A numpy array containing the indices of the nearest neighbors.
        """
        return nearest_neighbors_idx_njit(self._mat, vec, k)[0]

    def universal_markers(self, threshold: float = 0.95) -> npt.NDArray[np.uint32]:
        """
        Calculate (possible) universal markers based on the given count matrix.

        :param threshold: Fraction of reference sequences where a marker has to occur exatly once in order to be
        considered universal

        :return: The (column) indices of possible universal markers
        """
        return np.where(np.count_nonzero(self._mat == 1, axis=0) / self._mat.shape[0] >= threshold)[0]

    def estimate(
            self,
            vec: npt.NDArray[np.uint8],
            k: int,
            frac_eq: float = 1.0,
            knn_inds: Optional[npt.NDArray[np.uint64]] = None
    ) -> Tuple[float, float, int]:
        """
        Calculate an estimate for completeness and contamination for a vector based on common markers in the k nearest
        neighbors.
        :param vec: Input vector
        :param k: Number of nearest neighbors to consider. This parameter is ignored if knn_inds is provided.
        :param frac_eq: Fraction of similar counts within the nearest neighbors required to consider a Pfam as a
        marker
        :param knn_inds: If provided, this is used as the nearest neighbors (expects a 1D numpy array with database
        indices). If this is None, the function will determine the nearest neighbors itself.
        :return: A 3-tuple: First element is the completeness estimate, the second is the contamiation estimate
        (both between 0 and 1) and the third one is the number of markers that were used. (The last value is mainly
        intended for evaluation purposes.)
        """
        return estimate_njit(self.mat(), vec, k, frac_eq, knn_inds)


class QueryMatrix(Matrix[npt.NDArray[np.uint8]]):
    """
    A two-dimensional matrix containing Pfam counts. Contrary to the DatabaseMatrix this is intended to represent a set
    of CoCoPyE input bins/sequences.
    """
    _db_mat: Optional[npt.NDArray[np.uint8]] = None
    _db_metadata: Optional[pd.DataFrame] = None
    _k: Optional[int] = None
    _knn_inds: Optional[npt.NDArray[np.uint64]] = None
    _knn_scores: Optional[npt.NDArray[np.float32]] = None

    def __init__(self, mat: npt.NDArray[np.uint8]):
        """
        :param mat: 2-dimensional numpy matrix. Each row contains the counts of an input bin/sequence while each column
        represents a Pfam.
        """
        super().__init__(mat)

    def with_database(self, db: DatabaseMatrix, k: Optional[int] = None) -> QueryMatrix:
        """
        Add a DatabaseMatrix. This is required for almost all other functions of this class.

        :param db: DatabaseMatrix to be added.
        :param k: If provided, the function calculates and stores the k nearest neighbors for each query sequence. This
        is required for functions like estimates or taxonomy.
        """
        self._db_mat = db.mat()
        self._db_metadata = db.metadata()

        if k is not None:
            self._k = k
            self._knn_inds, self._knn_scores = nearest_neighbors_idx_njit_mat(self._db_mat, self._mat, self._k)

        return self

    def knn(self) -> Optional[Tuple[int, npt.NDArray[np.uint64]]]:
        """
        Return the previously calculated nearest neighbors. This requires a database matrix with k not None.

        :return: A Tuple containing the number of nearest neighbors as well as their indices. Returns None if this
        QueryMatrix doesn't contain any information about nearest neighbors.
        """
        if self._k is None:
            return None

        return self._k, self._knn_inds

    def knn_scores(self) -> Optional[npt.NDArray[np.float32]]:
        """
        This requires a database matrix with k not None.

        :return: Similarity scores for each query bin to its nearest neighbors (value between 0 and 1; a higher value
        indicates a higher similarity). Each row represents a bin, each column a neighbor.
        """
        if self._k is None:
            return None

        return self._knn_scores

    def preestimates(self, markers: npt.NDArray[np.uint32]) -> npt.NDArray[np.float32]:
        """
        Calculate preestimates based on universal markers.

        :param markers: An array containing (column) indices that will be used as markers

        :return: A two-dimensional numpy array. Each row contains two values which are (in this order) the completeness
        and the contamination estimate for the respective query bin.
        """
        submatrix = self._mat[:, markers]

        completeness = np.sum(np.clip(submatrix, 0, 1), axis=1) / markers.shape[0]
        contamination = np.sum(submatrix - np.clip(submatrix, 0, 1), axis=1) / markers.shape[0]

        return np.array([completeness, contamination], dtype=np.float32).T

    def estimates(
            self,
            frac_eq: float = 1.0,
            print_progress: bool = True
    ) -> Optional[npt.NDArray[np.float32]]:
        """
        Calculate a completeness and contamination estimate for all rows in the QueryMatrix based on common markers in
        the k nearest neighbors.

        This requires a database matrix with k not None.

        :param frac_eq: Fraction of similar counts within the nearest neighbors required to consider a Pfam as a marker
        :param print_progress: Print a progress bar to stdout.

        :return: A 2-dimensional numpy array. For each row in the QueryMatrix there is a row with two floats, where the
        first element ist the completeness and the second one the contamination estimate.
        """
        if self._db_mat is None:
            return None

        with ProgressBar(
                total=self.mat().shape[0],
                ncols=0,
                dynamic_ncols=False,
                desc="\033[0;37m[" + str(datetime.now()) + "]\033[0m Calculating estimates",
                disable=not print_progress
        ) as progress_bar:
            result = estimates_njit(self.mat(), self._db_mat, self._k, frac_eq, progress_bar, self._knn_inds)
        return result

    def taxonomy(self) -> Optional[List[Tuple[str, str]]]:
        """
        Use a consensus of the nearest neighbors for a taxonomy estimate. This will return the most specific taxonomic
        rank that is the same for all neighbors.

        :return: A list of strings representing taxonomy estimates, one for each query bin.
        """
        if self._db_metadata is None:
            return None

        metadata = self._db_metadata
        results = []

        for knn in self._knn_inds:
            found = False
            knn_meta = metadata.iloc[knn]
            for col in ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]:
                if knn_meta[col].isna().any():
                    continue
                if np.unique(knn_meta[col].to_numpy()).shape[0] == 1:
                    results.append((knn_meta.iloc[0][col], col))
                    found = True
                    break
            if not found:
                results.append(("-", "-"))

        return results

    def into_feature_mat(
            self,
            estimates: npt.NDArray[np.float32],
            resolution: int = 10
    ) -> Optional[FeatureMatrix]:
        if self._knn_inds is None:
            return None

        hist = Histogram(resolution)

        feature_vecs = []
        for idx, row in enumerate(self._mat):
            knn_mat = self._db_mat[self._knn_inds[idx]]

            knn_hists = []
            for neighbor in knn_mat:
                knn_hists.append(hist.calc_bins_for_two_counts(row, neighbor))

            mean_hist = np.mean(np.array(knn_hists), axis=0)
            mean_hist = mean_hist / np.sum(mean_hist)

            feature_vecs.append(np.concatenate([mean_hist, estimates[idx, :2]]))

        return FeatureMatrix(np.array(feature_vecs))


class FeatureMatrix(Matrix[npt.NDArray[np.double]]):
    """
    This class contains the features that are extracted from a QueryMatrix. They can be used to get completeness and
    contamination estimates of some machine learning model (in our case most likely a neural network).
    """
    def __init__(self, mat: npt.NDArray[np.double]):
        """
        You probably don't need to call this constructor directly. Instead you might want to call
        QueryMatrix.into_feature_mat.

        :param mat: Numpy matrix containing features
        """
        super().__init__(mat)

    def ml_estimates(self, model_file: str) -> npt.NDArray[np.float32]:
        """
        Calculate completeness or contamination estimates based on some machine learning model. (If you want both
        completeness and contamination, you have to call this function twice.)

        :param model_file: An sklearn .pickle file contaning the model that should be used.

        :return: An array contaning the predictions, one value for each query bin.
        """
        mlp = pickle.load(open(model_file, "rb"))
        return mlp.predict(self._mat)


def load_u8mat_from_file(filename: str) -> npt.NDArray[np.uint8]:
    """
    This is just a convenience function to load a numpy matrix from a file. If it doesn't suit your requirements, just
    use numpy.load or numpy.loadtxt directly.

    :param filename: Filename of the matrix file. If the extension is .npy it is assumed that the content is in binary
    format. If it is .npz the file will be treated as compressed binary format contaning exactly one matrix. Otherwise
    the file will be read as csv.

    :return: The loaded matrix as a numpy array
    """
    file_format = filename.split(".")[-1]

    mat: npt.NDArray[np.uint8]
    if file_format == "npy":
        mat = np.load(filename)
    elif file_format == "npz":
        mat = np.load(filename)["arr_0"]
    else:
        mat = np.loadtxt(filename, delimiter=",", dtype=np.uint8)

    return mat
