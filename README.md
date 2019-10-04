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

## Development

Clone higlass-python and higlass

```bash
git clone https://github.com/higlass/higlass
git clone https://github.com/higlass/higlass-python
```

Go to `higlass-python/js` and change the HiGlass dependency to `"higlass": "file:/absolute/path/to/higlass"` where `/absolute/path/to/higlass` should point to the directory you clone higlass into.

Finally, go to higlass-multicontact and activate its conda environment to do

```
jupyter labextension install ../path/to/higlass-python/js
```

After you have changed code in higlass you need to separately build higlass using `npm run build`.

After you have change code in higlass-jupyter you have to reinstall the labextension.
