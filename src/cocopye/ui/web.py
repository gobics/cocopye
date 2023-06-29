import os
import tempfile

import bottle
from bottle import route, run, request, view

from ..matrices import DatabaseMatrix, load_u8mat_from_file
from ..preprocessing.pfam import count_pfams


@route('/')
@view("upload")
def index():
    return dict()


@route('/upload', method='POST')
@view("upload")
def upload_file():
    fasta_file = request.files.get('fastaFile')

    name, ext = os.path.splitext(fasta_file.filename)
    if ext not in (".fasta", ".fna"):
        return dict(error="Unsupported file extension.")

    with tempfile.TemporaryDirectory() as tmpdirname:
        fasta_file.save(tmpdirname)

        config = request.app.config["config"]

        db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], "mat1234.npy")))
        query_mat, bin_ids = count_pfams(
            config["external"]["uproc_orf_bin"],
            config["external"]["uproc_bin"],
            config["external"]["uproc_db"],
            config["external"]["uproc_models"],
            tmpdirname,
            ext[1:]
        )

        estimates = query_mat.estimates(db_mat, 5)[0]

    return dict(result=(f'{estimates[0]*100:.2f}', f'{estimates[1]*100:.2f}'))


def start_server(config, debug=False):
    app = bottle.default_app()
    app.config["config"] = config

    if debug:
        run(host='localhost', port=8080, debug=True)
    else:
        run(host='localhost', port=8080, server="waitress", debug=False)
