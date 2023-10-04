# === HYPERPARAMETERS ==================================================================================================

K = 4
"""Number of nearest neighbors"""

RESOLUTION_COMP = 8
RESOLUTION_CONT = 12
"""Histogram resolution"""

# === DOWNLOAD URLS ====================================================================================================

UPROC = {
    "SRC": "http://uproc.gobics.de/downloads/uproc/uproc-1.2.0.tar.gz",
    "WIN": "http://uproc.gobics.de/downloads/uproc/uproc-1.2.0-win-x86_64.zip"
}
"""UProC source/binary download path"""

UPROC_MODEL = "http://uproc.gobics.de/downloads/models/model.tar.gz"
"""UProC model download path"""

PFAM_DB = {
    24: "http://uproc.gobics.de/downloads/db/pfam24.uprocdb.gz",
    28: "http://uproc.gobics.de/downloads/db/pfam28.uprocdb.gz"
}
"""UProC Pfam database download path"""

COCOPYE_DB = "https://github.com/gobics/cocopye-database/releases/latest/download/database.zip"
"""CocoPyE database download path"""

COCOPYE_DB_LATEST_RELEASE = "https://api.github.com/repos/gobics/cocopye-database/releases/latest"
