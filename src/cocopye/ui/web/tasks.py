import os

from celery import Celery

from ...matrices import DatabaseMatrix, load_u8mat_from_file
from ...preprocessing.pfam import count_pfams


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
def estimate_task(self, config, infolder: str):
    self.update_state(state="RUNNING")

    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], "mat1234.npy")))

    query_mat, bin_ids = count_pfams(
        config["external"]["uproc_orf_bin"],
        config["external"]["uproc_bin"],
        config["external"]["uproc_db"],
        config["external"]["uproc_models"],
        infolder,
    )

    return list(query_mat.estimates(db_mat, 5)[0])
