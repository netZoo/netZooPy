==========
Changelog
==========


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
