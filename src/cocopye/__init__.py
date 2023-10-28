"""
This is the API documentation of CoCoPyE.

### Submodules

- **constants**: Some constants that are required during runtime
- **core**: High-level functions that are used by `cocopye run`. If you are looking for a way to use the tool in a
Python script (without spawning a subprocess), this might be the module that you want to use.
- **histogram**: A histogram class that is used to convert a QueryMatrix into a FeatureMatrix
- **matrices**: Several matrix classes that are required during different stages of the tool
- **pfam**: Read FASTA files and use UProC to get their Pfam counts
- **ui**: Everything that is used for interaction with the user. Contains the command line interface, the webserver, the
automatic download of dependencies and configuration-related functions.
"""
