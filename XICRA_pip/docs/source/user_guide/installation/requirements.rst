.. ########################
.. _python-modules-required:
.. ########################

Python modules
--------------

There are several extra python module requirements that 
are needed and are summarized in the following table: 

.. csv-table::
   :header: "Module", "Version"
   :file: ../../../../XICRA/config/python/python_requirements.csv

These modules might have extra dependencies. Details of the list of 
all modules required are listed in :file:`XICRA/config/python/python_requirements.txt`. 
And accessible here:

.. toctree:: 
   python-requirements.rst

Although these dependencies will be fulfilled during the ``XICRA``
installation with ``pip``, you might be interested in installing them yourself. 

Using ``pip`` we can install them all at a glance. 

.. code-block:: sh

   pip install -r ./XICRA/config/python/python_requirements.txt

But again, following installation recommendations, we encourage you to create and install them 
within a virtual environment (See section: :ref:`Python environment<virtual-env-XICRA>` 
section for details).

You can test the presence of these ``python`` modules using the ``XICRA config`` module. 
Once you identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``PYTHONPATH`` or environment or you can execute the ``XICRA config`` 
with ``install`` option.

.. ######################
.. _soft-dependencies:
.. ######################

Software dependencies
---------------------

Also, several software packages are also required. They are listed in
:file:`XICRA/config/software/dependencies.csv`, which is shown below:

.. csv-table::
   :header-rows: 1 
   :file: ../../../../XICRA/config/software/soft_dependencies.csv

Most of the software are common software that any person doing bioinformatics should have, so
you might have already available within your system. However, installing ``XICRA`` using 
``conda`` all the dependecies will be correctky installed.   

You can test for any missing software dependencies using the ``XICRA config`` module. Once you 
identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``$PATH`` or you can execute the ``X config`` 
with ``install`` option.


.. #### Include links
.. include:: ../../links.inc