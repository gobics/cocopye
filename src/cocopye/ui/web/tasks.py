import os

from celery import Celery

from ...matrices import DatabaseMatrix, load_u8mat_from_file, QueryMatrix
from ...pfam import count_pfams


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
def estimate_task(self, config, pfam_version, infolder: str):
    self.update_state(state="RUNNING")

    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], pfam_version, "count_matrix.npz")))

    query_mat, bin_ids, _ = count_pfams(
        config["external"]["uproc_orf_bin"],
        config["external"]["uproc_prot_bin"],
        os.path.join(config["external"]["uproc_pfam_db"], pfam_version),
        config["external"]["uproc_models"],
        infolder,
    )

    query_mat = QueryMatrix(query_mat).with_database(db_mat, 4)

    return list(query_mat.estimates()[0])
