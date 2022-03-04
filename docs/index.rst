.. netZooPy documentation master file, created by
   sphinx-quickstart on Mon Oct  7 21:05:11 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

Home
===================================

netZooPy is a catalog of methods for the reconstruction and analysis of gene regulatory networks.

Available animals are: (gpu)PANDA, (gpu)LIONESS, (gpu)PUMA, SAMBAR, CONDOR, OTTER, and DRAGON.


netZooPy has merged of the following repositories and forks:

- `pypanda <https://github.com/davidvi/pypanda>`_ from David Van Ijzendorn, 

- Alessandro Marin's pypanda and Marieke Kuijjer's pypuma `fork <https://github.com/aless80/pypanda>`_, and 

- Cho-Yi Chen's `fork <https://github.com/QuackenbushLab/pypanda>`_.

- Genis Calderer's `pysambar <https://github.com/genis/pysambar>`_ and `pycondor <https://github.com/genis/pycondr>`_


Getting started
------------------

.. toctree::

   install/index

.. toctree::

   tutos/index

Documentation
-------------

Netzoopy also has a command line interface that allows to run the animal methods
from the command line, for example: 

.. code-block::

   netzoopy panda --e expression.txt --m motif.txt --p ppi.txt --o output_panda.txt

Check the documentation below for a full list of options.

.. toctree::
   functions/cli

.. toctree::

   functions/api


.. toctree::

   changelog



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
