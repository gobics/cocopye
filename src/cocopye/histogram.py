# Copyright 2023 Peter Meinicke, Niklas Birth, Nicolina Leppich
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

import numpy as np
import numpy.typing as npt

_MAX_COUNT = 255


class Histogram:
    """
    A Histogram class that is used for the transformation of a `cocopye.matrices.QueryMatrix` into a
    `cocopye.matrices.FeatureMatrix`.
    """
    _indx_mat: npt.NDArray[np.int32]
    _n_histogram_bins: int

    def __init__(self, resolution: int):
        """
        :param resolution: Histogram resolution. A higher value leads to more histogram bins and therefore a larger
        feature vector.
        """
        self._prepare_histogram(resolution)

    def _prepare_histogram(self, resolution: int) -> None:
        # calcutale borders of bins
        edge_vec = _calc_cr_hist_edges(resolution)
        n_histogram_bins = len(edge_vec) + 1
        # calculate lookup table for ratios
        indx_mat = _calc_index_matrix(edge_vec, n_histogram_bins)

        self._indx_mat = indx_mat
        self._n_histogram_bins = n_histogram_bins

    def calc_bins_for_two_counts(
            self, in_vec: npt.NDArray[np.uint8],
            neighbor_vec: npt.NDArray[np.uint8]
    ) -> npt.NDArray[np.int64]:
        """
        Calculate the bins for two vectors (usually Pfam counts of a query and a neighbor).

        :param in_vec: first (query) vector
        :param neighbor_vec: second (neighbor) vector

        :return: Histogram bins as numpy array
        """
        # clip to maximum count (255 for uint8)!
        # use both count vectors as indices
        # indvec contains all the ratio possitions for both count vect that would be in edge_vec
        ind_vec = self._indx_mat[np.clip(in_vec, 0, _MAX_COUNT), np.clip(neighbor_vec, 0, _MAX_COUNT)]
        # histogram counts for all bins
        x_new_vec = np.bincount(ind_vec[ind_vec != -1], minlength=self._n_histogram_bins)
        return x_new_vec


def _calc_cr_hist_edges(resolution: int) -> npt.NDArray[np.float32]:
    m = resolution + 1  # max. count hyperparameter
    # list of possible count ratios <= 1
    cr_list = [i / j for i in range(1, m) for j in range(i, m)]
    # sorted and non-redundant, now all possible ratios that can occcur
    cr_vec = np.unique(cr_list)

    # build difference to see if two ratios are too close together to form own bins
    diff_vec = cr_vec[1:] - cr_vec[:-1]
    # bulid edges of bins half way between the ratios
    sum_vec = (cr_vec[1:] + cr_vec[:-1]) / 2
    # avoid spurious bins due too numeric "errors"
    edge_half_vec = sum_vec[diff_vec > 1e-8]
    # np.histogram with edges does not provide open intervals for the left/right-most bins
    # we have to add extreme edges for artificial boundaries ...
    # [np.finfo('d').max] ==> Machine limits for floating point types
    # flip takes care of ratios above 1
    edge_vec = np.concatenate(
        (edge_half_vec, np.flip(1 / edge_half_vec)))

    return edge_vec


def _calc_index_matrix(edge_vec: npt.NDArray[np.float32], n_histogram_bins: int) -> npt.NDArray[np.int32]:
    # using 2D index table (256 x 256)!
    indx_mat = np.zeros((_MAX_COUNT + 1, _MAX_COUNT + 1), dtype=np.int32)
    # contains all ratios of 1-255 to 1-255
    # 1/1,  1/2,    1/3     ....    1/255
    # 2/1,  2/2,    2/3     ....    2/255
    # ...
    # 255/1,255/2,  255/3   ....  255/255
    cr_list = [i / j for i in range(1, _MAX_COUNT + 1)
               for j in range(1, _MAX_COUNT + 1)]
    # list of i and j indicies that belong to the same positions in ratios in cr_list
    cr_inds = [(i, j) for i in range(1, _MAX_COUNT + 1)
               for j in range(1, _MAX_COUNT + 1)]
    # tell us where in edge_vec to entries from cr_list need to go
    ind_vec = np.searchsorted(edge_vec, cr_list)
    for ind, val in enumerate(cr_list):
        i, j = cr_inds[ind]
        # tells us vor each count pair i and j with bin belongs to them
        indx_mat[i, j] = ind_vec[ind]

    indx_mat[:, 0] = n_histogram_bins - 1  # last histogram bin (index)
    # 0/0 is not an admissible count ratio -> provoke index error!
    indx_mat[0, 0] = -1
    return indx_mat
