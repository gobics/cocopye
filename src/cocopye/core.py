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

import os
from datetime import datetime
from typing import List, Dict

import numpy as np
import pandas as pd

from . import constants
from .matrices import DatabaseMatrix, load_u8mat_from_file, QueryMatrix
from .pfam import count_pfams


class Result:
    bin_id: str
    stage: int
    method: str
    count_ratio: float
    knn_scores: float
    taxonomy: str
    taxonomy_level: str
    notes: str

    comp_1_bac: float
    comp_1_arc: float
    cont_1_bac: float
    cont_1_arc: float

    comp_2: float
    cont_2: float
    num_markers_2: int

    comp_3: float
    cont_3: float

    def completeness(self):
        if self.stage == 1:
            return -1
        elif self.stage == 2:
            return self.comp_2
        else:
            return self.comp_3

    def contamination(self) -> float:
        if self.stage == 1:
            return -1
        elif self.stage == 2:
            return self.cont_2
        else:
            return self.cont_3

    def to_csv(self, verbosity: str = "standard") -> str:
        output_lists = {
            "standard": [self.bin_id, "{:.4f}".format(self.completeness()), "{:.4f}".format(self.contamination()), self.method, self.taxonomy,
                         self.taxonomy_level, self.notes],
            "extended": [self.bin_id, "{:.4f}".format(self.completeness()), "{:.4f}".format(self.contamination()), self.stage, self.method, self.num_markers_2,
                         self.count_ratio, self.knn_scores, self.taxonomy, self.taxonomy_level, self.notes],
            "full": [self.bin_id, self.stage, self.method, self.comp_1_arc, self.cont_1_arc, self.comp_1_bac, self.cont_1_bac,
                           self.comp_2, self.cont_2, self.num_markers_2, self.comp_3, self.cont_3, self.count_ratio,
                           self.knn_scores, self.taxonomy, self.taxonomy_level, self.notes]
        }

        return ",".join([str(item) for item in output_lists[verbosity]])

    def to_web(self) -> Dict[str, str]:
        return {
            "completeness": "{:.2f}%".format(self.completeness() * 100),
            "contamination": "{:.2f}%".format(self.contamination() * 100),
            "stage": self.stage,
            "taxonomy": self.taxonomy + " (" + self.taxonomy_level + ")",
        }


def core(cocopye_db: str,
         uproc_orf: str,
         uproc_prot: str,
         pfam_db: str,
         uproc_model: str,
         infolder: str,
         pfam_version: int,
         file_extension: str,
         num_threads: int,
         print_progress: bool = True
         ) -> List[Result]:
    pfam_version = str(pfam_version)

    log("Loading CoCoPyE database", print_progress)
    db_mat = DatabaseMatrix(
        load_u8mat_from_file(os.path.join(cocopye_db, pfam_version, "count_matrix.npz")),
        pd.read_csv(os.path.join(cocopye_db, pfam_version, "metadata.csv"), sep=",")
    )

    query_mat, bin_ids, count_ratio = count_pfams(
        uproc_orf,
        uproc_prot,
        os.path.join(pfam_db, pfam_version),
        uproc_model,
        infolder,
        file_extension,
        num_threads,
        print_progress
    )

    log("Determining nearest neighbors", print_progress)
    query_mat = QueryMatrix(query_mat).with_database(db_mat, constants.K)

    assert len(bin_ids) == query_mat.mat().shape[0]

    log("Calculating preestimates", print_progress)
    universal_arc = np.load(
        os.path.join(cocopye_db, pfam_version, "universal_Archaea.npy"))
    universal_bac = np.load(
        os.path.join(cocopye_db, pfam_version, "universal_Bacteria.npy"))

    preestimates_arc = query_mat.preestimates(universal_arc)
    preestimates_bac = query_mat.preestimates(universal_bac)

    estimates = query_mat.estimates(print_progress=print_progress, frac_eq=constants.FRAC_EQ)

    log("Calculating ML estimates", print_progress)
    feature_mat_comp = query_mat.into_feature_mat(estimates, constants.RESOLUTION_COMP[int(pfam_version)])
    feature_mat_cont = query_mat.into_feature_mat(estimates, constants.RESOLUTION_CONT[int(pfam_version)])

    ml_estimates_comp = feature_mat_comp.ml_estimates(
        os.path.join(cocopye_db, pfam_version, "model_comp.pickle")).clip(0, 1)
    ml_estimates_cont = feature_mat_cont.ml_estimates(
        os.path.join(cocopye_db, pfam_version, "model_cont.pickle")).clip(0, 1000000)

    log("Processing results", print_progress)
    taxonomy = query_mat.taxonomy()
    knn_scores = query_mat.knn_scores()

    completeness = []
    contamination = []
    stage = []
    for idx in range(len(bin_ids)):
        stage.append(1)
        if preestimates_bac[idx, 0] > preestimates_arc[idx, 0]:
            completeness.append(preestimates_bac[idx, 0])
            contamination.append(preestimates_bac[idx, 1])
        else:
            completeness.append(preestimates_arc[idx, 0])
            contamination.append(preestimates_arc[idx, 1])

        if completeness[idx] < constants.TRANSITION_1_2_MIN_COMP:
            continue

        stage[idx] = 2
        completeness[idx] = estimates[idx, 0]
        contamination[idx] = estimates[idx, 1]

        if completeness[idx] < constants.TRANSITION_2_3_MIN_COMP or contamination[idx] > constants.TRANSITION_2_3_MIN_CONT:
            continue

        stage[idx] = 3
        completeness[idx] = ml_estimates_comp[idx]
        contamination[idx] = ml_estimates_cont[idx]

    notes = []
    for idx in range(len(bin_ids)):
        # Currently we do not have any additional notes. But at least we could add some if we want.
        notes.append("")

    results = []
    for idx in range(len(bin_ids)):
        result = Result()

        result.bin_id = bin_ids[idx]
        result.stage = stage[idx]
        result.method = ["rejected", "markers", "markers + neural network"][stage[idx] - 1]
        result.count_ratio = count_ratio[idx]
        result.knn_scores = knn_scores[idx]
        result.taxonomy = taxonomy[idx][0]
        result.taxonomy_level = taxonomy[idx][1]
        result.notes = notes[idx]

        result.comp_1_bac = preestimates_bac[idx, 0]
        result.comp_1_arc = preestimates_arc[idx, 0]
        result.cont_1_bac = preestimates_bac[idx, 1]
        result.cont_1_arc = preestimates_arc[idx, 1]

        result.comp_2 = estimates[idx, 0]
        result.cont_2 = estimates[idx, 1]
        result.num_markers_2 = estimates[idx, 2]

        result.comp_3 = ml_estimates_comp[idx]
        result.cont_3 = ml_estimates_cont[idx]

        results.append(result)

    return results


def log(message: str, show: bool = True):
    """
    Wrapper around the print function for logging purposes.
    """
    if show:
        print("\033[0;37m[" + str(datetime.now()) + "]\033[0m " + message)
