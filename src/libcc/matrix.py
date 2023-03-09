from __future__ import annotations
import numpy as np
import numpy.typing as npt


class CountMatrix:
    _mat: npt.NDArray[np.uint8]

    def __init__(self, mat: npt.NDArray[np.uint8]):
        assert mat.ndim == 2, "Matrix has to be 2-dimensional"

        self._mat = mat

    def nearest_neighbors(self, vec: npt.NDArray[np.uint8], k: int) -> CountMatrix:
        return CountMatrix(self._mat[self.nearest_neighbors_idx(vec, k), :])

    def nearest_neighbors_idx(self, vec: npt.NDArray[np.uint8], k: int) -> npt.NDArray[np.int64]:
        num_refs, num_count = self._mat.shape

        assert vec.ndim == 1, "Vector has to be 1-dimensional"
        assert vec.shape[0] == num_count, "Vector length must be equal to the number of columns of the matrix"

        eq_count = np.zeros(num_refs)
        for idx, row in enumerate(self._mat):
            eq_count[idx] = np.sum(np.logical_and(0 < vec < 255, row == vec))

        return np.flip(np.argsort(eq_count))[:k]

    def transform_features(self, vec: npt.NDArray[np.uint8], n: int = 100, minmax: float = 3.) -> FeatureMatrix:
        eps = 1e-8
        n_rows = self._mat.shape[0]
        transformed_mat = np.zeros((n_rows, n))

        for i in range(n_rows):
            # Normal values (where at least one part is > 0 and we have no 255)
            inds = np.where(np.logical_and(np.logical_or(vec > 0, self._mat[i] > 0), vec < 255 and self._mat[i] < 255))[0]
            normal_values = np.log2((vec[inds] + eps) / (self._mat[i, inds] + eps))

            # Special values where at least one part contains max value (255)
            inds = np.where(np.logical_or(vec == 255, self._mat[i] == 255))[0]
            max_values = np.log2((vec[inds] + eps) / (self._mat[i, inds] + eps))

            # Count all the values and store them in our matrix
            transformed_mat[i, 0] += np.sum(np.concatenate((normal_values, max_values)) > minmax)
            transformed_mat[i, 1:] = np.histogram(normal_values, n - 2, range=(-minmax, minmax))[0]
            transformed_mat[i, -1] += np.sum(np.concatenate((normal_values, max_values)) < -minmax)

            # Do we need this? TODO
            # transformed_mat[i] /= np.sum(transformed_mat[i])

        return FeatureMatrix(transformed_mat)

    def save_to_file(self, filename: str) -> None:
        file_format = filename.split(".")[-1]

        if file_format == "npy":
            np.save(filename, self._mat)
        else:
            np.savetxt(filename, self._mat, delimiter=",")

    def __str__(self) -> str:
        return str(self._mat)


class FeatureMatrix:
    _mat: npt.NDArray[np.double]

    def __init__(self, mat: npt.NDArray[np.double]):
        assert mat.ndim == 2, "Matrix has to be 2-dimensional"

        self._mat = mat

    def mat(self) -> npt.NDArray[np.double]:
        return self._mat


def load_u8mat_from_file(filename: str) -> npt.NDArray[np.uint8]:
    file_format = filename.split(".")[-1]

    mat: npt.NDArray[np.uint8]
    if file_format == "npy":
        mat = np.load(filename)
    else:
        mat = np.loadtxt(filename, delimiter=",", dtype=np.uint8)

    return mat
