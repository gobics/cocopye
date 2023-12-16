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

import asyncio
import uuid
from datetime import datetime
from importlib import resources
import os
import subprocess
import werkzeug

from Bio import SeqIO
from fastapi import FastAPI, UploadFile, Request, WebSocket
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse
from jinja2 import Environment, PackageLoader, select_autoescape
import uvicorn
import shutil

from .. import config
from .tasks import estimate_task
from ..config import parse_config

CONFIG = parse_config()[1]

app = FastAPI()
jinja_env = Environment(
    loader=PackageLoader("cocopye.ui.web"),
    autoescape=select_autoescape()
)


@app.get("/", response_class=HTMLResponse)
async def root(request: Request):
    return jinja_env.get_template("index.html").render(
        request=request,
        contact=CONFIG["server"]["site_variables"]["contact"],
        limit=CONFIG["server"]["site_variables"]["upload_limit"],
        imprint=CONFIG["server"]["site_variables"]["imprint"],
        privacy=CONFIG["server"]["site_variables"]["privacy_policy"]
    )


@app.post("/upload")
async def upload(file: UploadFile):
    ws_id = str(uuid.uuid4())
    tmpdir = CONFIG["server"]["tmpdir"]

    os.makedirs(os.path.join(tmpdir, ws_id))

    with open(os.path.join(tmpdir, ws_id, werkzeug.utils.secure_filename(file.filename) + ".fna"), "wb") as outfile:
        # Not sure about this one. Maybe we should try to use some kind of non-blocking copy.
        shutil.copyfileobj(file.file, outfile)

    return {"ws_id": ws_id}


@app.websocket("/ws/{client_id}")
async def ws_endpoint(ws: WebSocket, client_id: str):
    await ws.accept()

    infolder = os.path.join(CONFIG["server"]["tmpdir"], client_id)
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

    print("Running webserver on http://" + config.CONFIG["server"]["host"] + ":" + str(config.CONFIG["server"]["port"]))

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
                "--without-heartbeat", "--without-gossip", "--without-mingle"], env=celery_env, stdout=subprocess.DEVNULL
        )

        uvicorn.run(
            "cocopye.ui.web.server:app",
            host=config.CONFIG["server"]["host"],
            port=config.CONFIG["server"]["port"],
            workers=config.CONFIG["server"]["workers"],
            log_level="error"
        )
