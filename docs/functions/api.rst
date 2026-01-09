The Animals (API)
-------------------

Animals listed in alphabetical order.

BONOBO
======

To run BONOBO you can first initialize the class with the expression data, 
and then run the compute_bonobo function to get the results: 

.. code-block:: 

   bonobo_obj_sparse = Bonobo(expression_file)
   bonobo_obj_sparse.run_bonobo(keep_in_memory=True, output_fmt='.hdf', sparsify=True, output_folder='../data/processed/bonobo_sparse_pvals/', save_pvals=False)

Here are the main functions and classes in the BONOBO module:

.. autoclass:: netZooPy.ligress.bonobo.Bonobo
   :members:
   :undoc-members:
   :inherited-members:
   :private-members:

Other functions
~~~~~~~~~~~~~~~

.. autofunction:: netZooPy.ligress.bonobo.compute_bonobo

.. autofunction:: netZooPy.ligress.io.prepare_expression

CONDOR
======

.. autoclass:: netZooPy.condor.condor.condor_object
   :members:
   :undoc-members:
   :inherited-members:
   :private-members:

Other function
~~~~~~~~~~~~~~~
.. autofunction:: netZooPy.condor.run_condor

LIONESS
=========


.. autoclass:: netZooPy.lioness.lioness.Lioness
   :members:
   :undoc-members:
   :inherited-members:
   :private-members:

.. autoclass:: netZooPy.lioness.analyze_lioness.AnalyzeLioness
   :members:
   :undoc-members:
   :private-members:





PANDA
======

.. autoclass:: netZooPy.panda.panda.Panda
   :members:
   :undoc-members:
   :inherited-members:
   :private-members:

.. autoclass:: netZooPy.panda.analyze_panda.AnalyzePanda
   :members:
   :undoc-members:
   :private-members:


PUMA
=====

.. autoclass:: netZooPy.puma.puma.Puma
   :members:
   :undoc-members:
   :inherited-members:
   :private-members:

SAMBAR
======

.. autofunction:: netZooPy.sambar.sambar

Extra sambar functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: netZooPy.sambar.corgenelength

.. autofunction:: netZooPy.sambar.convertgmt

.. autofunction:: netZooPy.sambar.desparsify

.. autofunction:: netZooPy.sambar.binomial_dist

.. autofunction:: netZooPy.sambar.clustering
