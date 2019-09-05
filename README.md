# HiGlass Multicontact

Notebooks for experimenting with multicontact data.

## Getting started

#### Requirements

- Conda `v4.6.14` (or higher)
- FUSE: [libfuse](https://github.com/libfuse/libfuse) `v2.9.7` (or higher) or [OSXFUSE](https://osxfuse.github.io/) `v3.9.2` (or higher)

_Other versions might work too but we only tested the above mentioned versions._

#### Installation

First, install the environment:

```
git clone https://github.com/4dn-dcic/higlass-multicontact
cd higlass-multicontact
conda env create -f environment.yml
```

Next, install HiGlass' jupyter extension:

```
conda activate higlass-multicontact
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install higlass-jupyter
```

Finally, start Jupyterlab:

```
jupyter-lab
```
