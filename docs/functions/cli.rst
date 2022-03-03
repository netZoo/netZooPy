Command Line Interface
-----------------------


Netzoopy also has a command line interface that allows to run the animal methods
from the command line, for example: 

.. code-block::

   netzoopy panda --e expression.txt --m motif.txt --p ppi.txt --o output_panda.txt

The command above would run the panda message passing method using the expression, 
motif and ppi data passed as input, and returning an output file containing the 
panda weighted edges.

While we expect to provide most animals from CLI, so far you can use: panda,lioness.

If you want to use the previous command line calls (python run_panda.py), check the 
:ref:`Run functions (legacy)`

Commands
=========================

.. click:: netZooPy.cli:cli
   :prog: netzoopy
   :commands: panda,lioness,condor
   :nested: full


Run functions (legacy)
=========================

Panda
~~~~~~~

.. autofunction:: netZooPy.panda.run_panda.main

Lioness
~~~~~~~~~

.. autofunction:: netZooPy.lioness.run_lioness.main
