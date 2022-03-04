# Installation guide

You can install netZooPy through pip or through conda. Below, you will find 
all the steps and requirements for both cases.

## PIP installation

### Requirements

- Python 3

In addition to the following pip packages:

- networkx

- numpy

- pandas

- matplotlib

- scipy

- python-igraph

### Install

- `git clone https://github.com/netZoo/netZooPy.git`

- `cd netZooPy`

- `pip3 install -e .`

- Then you can import netZooPy in your code through `import netZooPy`

### Troubleshooting

- If you experience installation issues with igraph on MacOs. Please install igraph through `brew install igraph`.

  Please refer to the detailed error fix for MacOs and Ubuntu [here](https://github.com/igraph/python-igraph/issues/159).

  Refer to the igraph instrallation [manual](https://igraph.org/python/#pyinstallosx)

- To report any installation issue or function bug, please report through opening an [issue](https://github.com/netZoo/netZooPy/issues) on github.


## Conda installation

On anaconda.org you will find the conda recipes for all platforms. We recommend using conda environments to keep your analyses self-contained and reproducible.

To install netzoopy through conda:

```bash
conda install -c netzoo -c conda-forge netzoopy
```