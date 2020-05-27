# Branch Description

Given the inability for R to use script `panda.py` and script `lioness.py` located in master branch [here](https://github.com/netZoo/netZooPy/blob/master/netZooPy/panda/panda.py) and [here](https://github.com/netZoo/netZooPy/blob/master/netZooPy/lioness/lioness.py) casued by importing module error (i.e. R package reticulate can not process importing modules from other sciprts, like `from .timer import Timer` and `from netZooPy.panda.panda import Panda` will cause an erorr `ValueError: Attempted relative import in non-package`), a specific branch created here is used in the R package [netZooR](https://github.com/netZoo/netZooR) to source the Python implementation of PANDA and LIONESS.

## Details

- `panda.py` script here derived from class `Panda` and class `Timer` in the Python library [`netZooPy`](https://github.com/netZoo/netZooPy)(i.e. the master branch of this repository) is used in function [panda.py()](https://github.com/netZoo/netZooR/blob/master/R/PANDA.R) in the netZooR package.

- `lioness.py` script here derived from class `Lioness` and class `Timer` in Python library `netZooPy` is used in function [lioness.py()](https://github.com/netZoo/netZooR/blob/master/R/LIONESS.R) in the netZooR package.




## Development reminder

Please keep the content of the scirpts in this branch consistent with up-to-date version of netZooPy. 

Current contenet is based on PR[#96](https://github.com/netZoo/netZooPy/pull/96)
