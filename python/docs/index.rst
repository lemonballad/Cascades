Cascades Documentation
======================

Cascade artifact simulations for 2D resonance Raman spectroscopy.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   api

Installation
------------

.. code-block:: bash

   cd python
   pip install -e .

Quick Start
-----------

.. code-block:: python

   from cascades.simulations import run_2drr_simulation

   results = run_2drr_simulation(
       solvent="methanol",
       nmode=3,
       nquanta=4
   )

API Reference
-------------

Core Modules
~~~~~~~~~~~~

.. automodule:: cascades.core.basis
   :members:
   :undoc-members:

.. automodule:: cascades.core.franck_condon
   :members:
   :undoc-members:

.. automodule:: cascades.core.response
   :members:
   :undoc-members:

.. automodule:: cascades.core.fsrs
   :members:
   :undoc-members:

.. automodule:: cascades.core.offres
   :members:
   :undoc-members:

Parameters
~~~~~~~~~~

.. automodule:: cascades.parameters.pna
   :members:
   :undoc-members:

.. automodule:: cascades.parameters.myoglobin
   :members:
   :undoc-members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
