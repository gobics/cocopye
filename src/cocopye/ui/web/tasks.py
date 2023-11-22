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
from typing import Dict

from celery import Celery

from ... import core


app = Celery(
    "tasks",
    backend="redis://localhost",
    broker="redis://localhost",
    # This does not work (crashes on task execution) and I have no idea why.
    # backend=os.getenv("BACKEND"),
    # broker=os.getenv("BROKER"),
    broker_connection_retry=False,
    broker_connection_retry_on_startup=True,
    broker_connection_max_retries=10,
)


@app.task(bind=True, time_limit=os.getenv("CELERY_TIME_LIMIT"))
def estimate_task(self, config, pfam_version: int, infolder: str, debug: bool = False) -> Dict[str, str]:
    self.update_state(state="RUNNING")

    result = core.core(config["external"]["cocopye_db"],
                       config["external"]["uproc_orf_bin"],
                       config["external"]["uproc_prot_bin"],
                       config["external"]["uproc_pfam_db"],
                       config["external"]["uproc_models"],
                       infolder,
                       pfam_version,
                       "fna",
                       1,
                       print_progress=debug
                       )[0]

    return result.to_web()
