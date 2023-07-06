import os

from celery import Celery

from ...matrices import DatabaseMatrix, load_u8mat_from_file
from ...preprocessing.pfam import count_pfams


app = Celery("tasks", backend="redis://localhost", broker="redis://localhost")


@app.task(bind=True)
def estimate(self, config, infolder: str):
    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], "mat1234.npy")))

    self.update_state(state="RUNNING", meta={"stage": 0})
    query_mat, bin_ids = count_pfams(
        config["external"]["uproc_orf_bin"],
        config["external"]["uproc_bin"],
        config["external"]["uproc_db"],
        config["external"]["uproc_models"],
        infolder,
    )

    self.update_state(state="RUNNING", meta={"stage": 1})
    return list(query_mat.estimates(db_mat, 5)[0])
