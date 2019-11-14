# pyCONDOR
Python implementation of the BRIM algorithm for bipartite community structure detection as described in "Modularity and community detection in bipartite networks" by Michael J. Barber." This package is somewhat translated from the R package CONDOR https://github.com/jplatig/condor described in the paper "Bipartite Community Structure of eQTLs" by John Platig , Peter J. Castaldi, Dawn DeMeo, John Quackenbush.

## Install instructions
This package uses the python libraries numpy,pandas and python-igraph (igraph).

The easiest way to install pyCONDOR is using pip:
```
git clone https://github.com/genisott/pycondor.git
cd pycondor
pip install .
```
To use the package you can import it using: ```import condor```. And then access the different functions implemented with ```condor.function()```.


## Usage instructions
Supose you have a network (weighted or not) as an edgelist loaded into a pandas dataframe.
```
import condor
co = condor.condor_object(net)
```
Returns a condor object (dictionary with a graph, node names, community membership dataframes...)
```
co = condor.initial_community(co)
```
Computes the initial community structure and updates the condor object.
```
co = condor.brim(co)
```
Runs the iterative modularity optimization algorithm (BRIM) and updates the condor object with the final membership.
To see the results type:
```
co["tar_memb"]
co["reg_memb"]
```
To compute the qscores for the vertices type:
```
co = condor.qscores(co) # Computes the qscores
co["qscores"]["reg_qscores"] # Dataframe containing the qscores for the regulators.
co["qscores"]["tar_qscores"] # Dataframe containing the qscores for the targets.
```
