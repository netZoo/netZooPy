==========
Changelog
==========


0.9.12
-------------------------

- We are changing the PANDA outputs and default flags. For now we are updating the command line 
 call only, behavior is kept as in 0.9.11 for the internal functions. By passing `old_compatible = False`
 the final output will always have column headers and indices.
- PUMA and PANDA do not save_tmp as default.
- lioness for puma has been fixed
- Fixed PANDA data preprocessing bug

0.9.11 (2022-11-04)
-------------------------

- Added LIONESS for DRAGON with tests
- PANDA preprocessing expression: In Panda preprocessing there was a problem with indices. Using gene2idx.get(x, 0) always give you the index 0 if x is missing fro
m gene2idx.get (like a gene in gene expression and not in motif, since gene2idx is build on top of the intersection of expression and motif). Now we use gene_names to
both create the indices for self.expression and to access with .loc[] the expression data frame self.expression_data
- New PANDA tests
- Updated LIONESS start and end parameters so that they are independent of the background. Example: One can now run panda on 100 samples
  and then apply LIONESS on only the first 10.
- Added LIONESS subset parameter: passing subset parameters (a list of indices or sample names, [1,2,10]) allows to run
  LIONESS only on specific samples. This parameter has priority over the start and end parameters.
  
0.9.10 (2022-10-28)
------------------

- Fixing single/double precision for GPU
- Clearing GPU after computation to free more memory

0.9.9 (2022-10-21)
------------------

- added the case for square nonsymmetric matrices for normalization in panda
- Updated tests for panda and lioness to match MATLAB
- Fixed Panda-Lioness GPU inconsistencies
- Forcing igraph<0.10, otherwise community assignment results change. This will need further investigation for the future.
- Fixed lioness GPU export (now lioness allows to save the full matrix, with explicit edge and sample names).

0.9.6 (2022-06-10)
------------------

- Ligress filters PPI according to input motif


0.9.5 (2022-05-24)
------------------

- Added output with sample names in Lioness
- ligress sample names are setup as strings
- correct order of motif prior in ligress

0.9.4 (2022-05-20)
------------------

- First ligress release
- solved puma bug


0.9.2 (2022-03-04)
------------------
- added command line interface (panda, lioness)
- updating documentation

0.9.0 (2022-02-11)
------------------
- we fixed the panda-lioness and puma-lioness behavior ( panda was passing the updated motif to lioness ). The results are now compatible with the ones of netzooR.
- removed py3.6 support
- updated version on anaconda.org

0.8.0 (2021-06-08)
-------------------
- support for Python v3.9 
- addition of DRAGON + unit tests +tutorial and many bug fixes that Daniel and Marouen have been doing as a user requests

0.7.2 (2020-07-18)
------------------

- PANDA reads arguments as dataframes in addition to file paths
- changed condor ground truth to match output of `python-igraph 0.8.2 <https://github.com/netZoo/netZooPy/issues/82>`_. 

0.7.1 (2020-06-27)
------------------

- Major fix for OTTER behavior across platforms.

0.7.0 (2020-01-18)
------------------

- new tool: OTTER
- unit test for OTTER
- fix for PANDA `force` field
- tweaks for compatibility of gpuPANDA with cupy

0.6.2 (Stockholm) (2020-05-15)
------------------------------

- Added gpuPANDA, which is a gpu-accelerated implementation of PANDA
- Added gpuLIONESS
- Added a gpuPANDA and gpuLIONESS tutorial
- Fixed condor dependency to python-igraph (still under investigation in #82 )

0.6.1 (2020-01-18)
------------------

- sambar tutorial
- condor tutorial
- added 3.8 to Ubunutu test server (along with 3.6 and 3.7 )
- Created three options for data processing in PANDA.
     - Union: adds rows for genes/TFs that are missing in at least one prior (expression, ppi, motif)
     - Intersection: removes TF/genes that missing in at least one prior
     - Legacy: previous data processing behavior
- The default was set to union in netZooM, netZooR, netZooPy as it is the default in netZooC.

0.5.0 (2019-11-22)
------------------

- pysambar

0.4.0 (2019-11-18)
------------------

- pycondor

0.3.0 (2019-11-14)
------------------

- pypuma

0.2.0 (2019-11-13)
------------------

- pylioness

0.1.1 (2019-9-3)
------------------

- fixed call to save_memory=True

0.1.0 (2019-7-26)
------------------

- transition to python 3
- Changelog added to the doc
- pypanda: original import and NaN values in normalized matrices are replaced with values normalized by the overall z-score. This allows running the Toy Data provided in this repository.  
