# === HYPERPARAMETERS ==================================================================================================

K = 4
"""Number of nearest neighbors"""
FRAC_EQ = 1.0
"""Fraction of counts that have to be the same in all neighbors for a Pfam to be considered as a marker"""

RESOLUTION_COMP = {24: 8, 28: 6}
"""Histogram resolution for completeness"""
RESOLUTION_CONT = {24: 7, 28: 12}
"""Histogram resolution for contamination"""

TRANSITION_1_2_MIN_COMP = 0.2
"""Minimal completeness to allow a bin to move from stage 1 to 2"""
TRANSITION_2_3_MIN_COMP = 0.55
"""Minimal completeness to allow a bin to move from stage 2 to 3"""

# === DOWNLOAD URLS ====================================================================================================

UPROC = {
    "SRC": "http://uproc.gobics.de/downloads/uproc/uproc-1.2.0.tar.gz",
    "WIN": "http://uproc.gobics.de/downloads/uproc/uproc-1.2.0-win-x86_64.zip"
}
"""
UProC source/binary download path.
Please note that there are some hardcoded directory names in cocopye.ui.external.bin. You might need to adjust them
if you change this URLs (for example in the unlikely case that there will be a UProC update).
"""

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
"""Latest release of the CoCoPyE database (GitHub API). This is used for update checking."""
