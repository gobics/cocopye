import asyncio
from importlib import resources
import os
import subprocess

from Bio import SeqIO
from fastapi import FastAPI, UploadFile, WebSocket
from fastapi.staticfiles import StaticFiles
import uvicorn
import shutil
import random

from ..config import CONFIG, ARGS
from .tasks import estimate_task

app = FastAPI()


@app.post("/upload")
async def upload(file: UploadFile):
    ws_id = random.randint(0, 100000000)
    tmpdir = CONFIG["server"]["tmpdir"]

    os.makedirs(os.path.join(tmpdir, str(ws_id)))

    with open(os.path.join(tmpdir, str(ws_id), file.filename + ".fna"), "wb") as outfile:
        # Not sure about this one. Maybe we should try to use some kind of non-blocking copy.
        shutil.copyfileobj(file.file, outfile)

    return {"ws_id": ws_id}


@app.websocket("/ws/{client_id}")
async def ws_endpoint(ws: WebSocket, client_id: int):
    await ws.accept()

    infolder = os.path.join(CONFIG["server"]["tmpdir"], str(client_id))
    infile = os.listdir(infolder)[0]

    await ws.send_text("Checking input file")
    for record in SeqIO.parse(os.path.join(infolder, infile), "fasta"):
        if any(c not in "ACGTN" for c in record.seq):
            await ws.send_text("Error: Invalid input file.")
            await ws.close()
            shutil.rmtree(infolder)
            return

    task = estimate_task.delay(CONFIG, ARGS.pfam_version, infolder)

    if task.state == "PENDING":
        await ws.send_text("Waiting for task execution")

    while task.state == "PENDING":
        # This is the best solution I found. There does not seem to be somthing like await estimates and estimates.get()
        # is a blocking call which freezes the webserver until the task is completed. (I'm actually not sure why this is
        # the case. I would have thought that it should be possible to work on three requests simultaneously if we start
        # uvicorn with three workers, but alredy the second one freezes.)
        await asyncio.sleep(0.3)

    await ws.send_text("Running CoCoPyE")

    while task.state == "RUNNING":
        await asyncio.sleep(0.3)

    if task.state == "SUCCESS":
        estimates = task.get()
        await ws.send_text(
            "Completeness: " + f'{estimates[0] * 100:.2f}%' + ", "
            "Contamination: " + f'{estimates[1] * 100:.2f}%'
        )
    else:
        await ws.send_text("Something went wrong.")

    await ws.close()
    shutil.rmtree(infolder)


app.mount(
    "/",
    StaticFiles(directory=str(resources.files("cocopye.ui.web.static").joinpath("")), html=True), name="static"
)


def run_server():
    # It seems to be quite difficult to get configuration options into the tasks module when it is called by celery.
    # Environment variables are probably not the best solution, but at least they work (or they don't; TODO).
    celery_env = os.environ.copy()
    # celery_env["CELERY_BROKER_URL"] = CONFIG["server"]["celery"]["broker"]
    # celery_env["CELERY_RESULT_BACKEND"] = CONFIG["server"]["celery"]["backend"]
    celery_env["CELERY_TIME_LIMIT"] = str(CONFIG["server"]["celery"]["time_limit"])
    if celery_env["CELERY_TIME_LIMIT"] == "0":
        celery_env["CELERY_TIME_LIMIT"] = "10000000"

    if CONFIG["server"]["debug"]:
        subprocess.Popen(
            ["celery", "-A", "cocopye.ui.web.tasks", "worker", "-c", str(CONFIG["server"]["celery"]["workers"])],
            env=celery_env
        )

        uvicorn.run(
            "cocopye.ui.web.server:app",
            host=CONFIG["server"]["host"],
            port=CONFIG["server"]["port"],
            workers=CONFIG["server"]["workers"],
        )
    else:
        subprocess.Popen(
            ["celery", "-A", "cocopye.ui.web.tasks", "worker", "-c", str(CONFIG["server"]["celery"]["workers"]),
                "--without-heartbeat", "--without-gossip", "--without-mingle"], env=celery_env
        )

        uvicorn.run(
            "cocopye.ui.web.server:app",
            host=CONFIG["server"]["host"],
            port=CONFIG["server"]["port"],
            workers=CONFIG["server"]["workers"],
            log_level="error"
        )
