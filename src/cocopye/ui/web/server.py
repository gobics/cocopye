import asyncio
from datetime import datetime
from importlib import resources
import os
import subprocess

from Bio import SeqIO
from fastapi import FastAPI, UploadFile, WebSocket
from fastapi.staticfiles import StaticFiles
import uvicorn
import shutil
import random

from .. import config
from .tasks import estimate_task
from ..config import parse_config

CONFIG = parse_config()[1]

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

    await ws.send_json({"status": "progress", "content": "Checking input file"})
    for record in SeqIO.parse(os.path.join(infolder, infile), "fasta"):
        if any(c not in "ACGTN" for c in record.seq):
            await ws.send_json({"status": "error", "content": "Invalid input file."})
            await ws.close()
            shutil.rmtree(infolder)
            return

    log("Starting new task")
    task = estimate_task.delay(
        CONFIG,
        24 if os.environ["COCOPYE_PFAM24"] == "1" else 28,
        infolder,
        CONFIG["server"]["debug"]
    )

    if task.state == "PENDING":
        await ws.send_json({"status": "progress", "content": "Waiting for task execution"})

    while task.state == "PENDING":
        # This is the best solution I found. There does not seem to be somthing like await estimates and estimates.get()
        # is a blocking call which freezes the webserver until the task is completed. (I'm actually not sure why this is
        # the case. I would have thought that it should be possible to work on three requests simultaneously if we start
        # uvicorn with three workers, but alredy the second one freezes.)
        await asyncio.sleep(0.3)

    await ws.send_json({"status": "progress", "content": "Running CoCoPyE"})

    while task.state == "RUNNING":
        await asyncio.sleep(0.3)

    if task.state == "SUCCESS":
        result = task.get()
        await ws.send_json({"status": "result", "content": result})
        log("Successfully completed task")
    else:
        await ws.send_json({"status": "error", "content": "Something went wrong."})
        log("Failed task")

    await ws.close()
    shutil.rmtree(infolder)


app.mount(
    "/",
    StaticFiles(directory=str(resources.files("cocopye.ui.web.static").joinpath("")), html=True), name="static"
)


def log(message: str):
    f = open(os.path.join(CONFIG["server"]["logdir"], "server.log"), "a")
    f.write("[" + str(datetime.now()) + "] " + message + "\n")
    f.close()


def run_server():
    os.environ["COCOPYE_PFAM24"] = "1" if config.ARGS.pfam24 else "0"
    os.makedirs(CONFIG["server"]["logdir"], exist_ok=True)

    # It seems to be quite difficult to get configuration options into the tasks module when it is called by celery.
    # Environment variables are probably not the best solution, but at least they work.
    celery_env = os.environ.copy()
    # ...except when they don't. But this is something we can fix later as it is not that important.
    # celery_env["CELERY_BROKER_URL"] = CONFIG["server"]["celery"]["broker"]
    # celery_env["CELERY_RESULT_BACKEND"] = CONFIG["server"]["celery"]["backend"]
    celery_env["CELERY_TIME_LIMIT"] = str(config.CONFIG["server"]["celery"]["time_limit"])
    if celery_env["CELERY_TIME_LIMIT"] == "0":
        celery_env["CELERY_TIME_LIMIT"] = "10000000"

    if config.CONFIG["server"]["debug"]:
        subprocess.Popen(
            ["celery", "-A", "cocopye.ui.web.tasks", "worker", "-c", str(config.CONFIG["server"]["celery"]["workers"])],
            env=celery_env
        )

        uvicorn.run(
            "cocopye.ui.web.server:app",
            host=config.CONFIG["server"]["host"],
            port=config.CONFIG["server"]["port"],
            workers=config.CONFIG["server"]["workers"],
        )
    else:
        subprocess.Popen(
            ["celery", "-A", "cocopye.ui.web.tasks", "worker", "-c", str(config.CONFIG["server"]["celery"]["workers"]),
                "--without-heartbeat", "--without-gossip", "--without-mingle"], env=celery_env
        )

        uvicorn.run(
            "cocopye.ui.web.server:app",
            host=config.CONFIG["server"]["host"],
            port=config.CONFIG["server"]["port"],
            workers=config.CONFIG["server"]["workers"],
            log_level="error"
        )
