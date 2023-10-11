import os
from typing import Dict

from celery import Celery

from ... import core


app = Celery(
    "tasks",
    backend="redis://localhost",
    broker="redis://localhost",
    # backend=os.getenv("BACKEND"),     This does not work (crashes on task execution) and I have no idea why.
    # broker=os.getenv("BROKER"),       TODO
    broker_connection_retry=False,
    broker_connection_retry_on_startup=True,
    broker_connection_max_retries=10,
)


@app.task(bind=True, time_limit=os.getenv("CELERY_TIME_LIMIT"))
def estimate_task(self, config, pfam_version: int, infolder: str) -> Dict[str, str]:
    self.update_state(state="RUNNING")

    result = core.core(config["external"]["cocopye_db"],
                       config["external"]["uproc_orf_bin"],
                       config["external"]["uproc_prot_bin"],
                       config["external"]["uproc_pfam_db"],
                       config["external"]["uproc_models"],
                       infolder,
                       pfam_version,
                       "fna",
                       1
                       )[0]

    return result.to_web()
