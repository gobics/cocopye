# CoCoPyE
**Documentation**:  [https://n.birth.pages.gwdg.de/CoCoPyE](https://n.birth.pages.gwdg.de/CoCoPyE)

*You can ignore the certificate error that occurs when opening the documentation.
This an issue with GitLab Pages that is caused by the dot in my username. It won't be
a problem anymore once we move this project to GitHub.*

## Installation

As soon as we publish the tool, it should be installable via PyPI by just running
`pip install cocopye`. However, currently there are some more steps required.

```
# Clone this repository and change directory
git clone https://gitlab.gwdg.de/n.birth/CoCoPyE.git
cd CoCoPyE

# create venv
python -m venv .venv

# activate venv (use the command for your OS)
source .venv/bin/activate  # Linux
.venv/Scripts/activate     # Windows

# install CoCoPyE package
pip install .
```

## Usage

```
cocopye run -i <bin folder> -o <outfile>
```

On the first run, CoCoPyE needs to download some external files. This won't be necessary on subsequent runs.

(I will write a more detailed usage guide in the future. but for now this is the most importang command.
For everything else, you can use the help flag: `cocopye -h` / `cocopye run -h`.)