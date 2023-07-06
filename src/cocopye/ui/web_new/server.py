import os
import subprocess

from Bio import SeqIO
from fastapi import FastAPI, UploadFile, Request, WebSocket
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse
import uvicorn
import shutil
import random

from ..config import CONFIG
from .tasks import estimate

app = FastAPI()
templates = Jinja2Templates(directory="views")


@app.get("/", response_class=HTMLResponse)
async def root(request: Request):
    return templates.TemplateResponse("main.html", { "request": request })


@app.post("/upload")
async def upload(file: UploadFile):
    ws_id = random.randint(0, 100000000)

    os.mkdir(os.path.join("/home/nemo/ExcludeFromBackup/CoCoPyE/files", str(ws_id)))

    with open(os.path.join("/home/nemo/ExcludeFromBackup/CoCoPyE/files", str(ws_id), file.filename + ".fna"), "wb") as outfile:
        shutil.copyfileobj(file.file, outfile)

    return {"ws_id": ws_id}


@app.websocket("/ws/{client_id}")
async def ws_endpoint(ws: WebSocket, client_id: int):
    await ws.accept()

    infolder = os.path.join("/home/nemo/ExcludeFromBackup/CoCoPyE/files", str(client_id))
    infile = os.listdir(infolder)[0]

    await ws.send_text("Checking input file")
    if not check_fasta(os.path.join(infolder, infile)):
        await ws.send_text("Error: Invalid input file.")
        await ws.close()
        shutil.rmtree(infolder)
        return

    async def on_raw_message(body):
        await ws.send_text(body["status"])

    estimates = estimate.delay(CONFIG, infolder).get(on_message=on_raw_message)

    await ws.send_text("Completeness: " + f'{estimates[0]*100:.2f}%' + ", Contamination: " + f'{estimates[1]*100:.2f}%')

    await ws.close()
    shutil.rmtree(infolder)


def check_fasta(infile):
    for record in SeqIO.parse(infile, "fasta"):
        if any(c not in "ACGTN" for c in record.seq):
            return False

    return True


def run_server():
    subprocess.Popen(
        ["celery", "-A", "cocopye.ui.web_new.tasks", "worker"]
    )

    uvicorn.run("cocopye.ui.web_new.server:app", host='127.0.0.1', port=8000, workers=3)
