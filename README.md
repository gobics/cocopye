# CoCoPyE

**CoCoPyE is a fast tool for quality assessment of microbial genomes. It is able to reliably predict
completeness and contamination of bacterial and archaeal genomes. Additionally, it can provide a
taxonomic classification of the input.**

Background: The classical approach for estimation of quality indices solely relies on a relatively small
number of universal single copy genes. Because these classical markers only cover a small fraction of the
whole genome the quality assessment can be rather unreliable. Our method is based on a novel
two-stage feature extraction and transformation scheme. It first performs a flexible extraction
of genomic markers and then refines the marker-based estimates with a machine learning approach based on
count-ratio histograms. In our simulation studies CoCoPyE showed a more accurate prediction of  quality
indices than existing tools.

## Getting started

CoCoPyE is available via pip ~~and conda (conda-forge channel).~~ Coming soon. See the [project wiki](https://github.com/gobics/cocopye/wiki)
for installation and usage instructions.

- [Quickstart](https://github.com/gobics/cocopye/wiki/Quickstart)
- [Installation](https://github.com/gobics/cocopye/wiki/Installation)
- [Usage](https://github.com/gobics/cocopye/wiki/Usage)

### Online Demo

You can test CoCoPyE without installation on [our project homepage](https://cocopye.uni-goettingen.de). Please note that the online demo can process only 
one query genome per request and is less performant than a local installation. Therefore it is highly recommended to use the online
version only for evaluation purposes and install CoCoPyE on your own machine for productive use.

## Additional notes

### Contact

For bug reports, suggestions or questions, please open an issue on [GitHub](https://github.com/gobics/cocopye/issues)
or send an email to [email@example.com](mailto:email@example.com). ***(TODO: Replace by the actual address.)***

### API documentation

You can find the API documentation of the CoCoPyE package on [https://gobics.github.com/cocopye](https://gobics.github.com/cocopye).

### License

CoCoPyE is available under the terms of the GNU General Public License, version 3 or later. See [`COPYING`](https://github.com/gobics/cocopye/blob/master/COPYING) for the full license text.
