.. ########################
.. _installing:
.. ########################

Installation
************

.. toctree::
   :hidden:
   
   requirements.rst
   python-environment.rst
      
.. contents::

This is an installation guide for the ``XICRA`` pipeline. 

``XICRA`` is a pipeline composed of multiple lines of code 
that calls and automatizes the analysis from multiple bioinformatic 
tools. Several dependencies (python, perl, additional software and their dependencies 
and third-party software) are required for the ``XICRA`` analysis. 

You first need to get ``XICRA`` (from different sources available): 
you can either install ``XICRA`` latest version using Conda_ and `Python pip`_ 
or if you want to contribute or see additional details of the project, it's 
recommended you :ref:`install the latest development version<install-from-source>`.

Once you get the code, before running ``XICRA`` you must make sure 
you have all the dependencies fulfilled from section 
:ref:`Requirements and dependencies<Requirements-dependencies>` using the 
``XICRA config`` module.


.. ########################
.. _Requirements-dependencies:
.. ########################

Requirements and dependencies
=============================

.. ########################
.. _system-requirements:
.. ########################

System requirements
-------------------

Take into account that you may need some system requirements to install ``XICRA`` already available
in your system such as ``python3``, ``python3-dev`` and ``python3-venv`` or ``build-essential`` libraries among others. 

Check and install them by typing:

.. code-block:: sh
   
   sudo apt install python3 python3-dev python-dev python3-venv python-pip
   sudo apt install build-essential libssl-dev libffi-dev libxml2-dev libxslt1-dev zlib1g-dev

``XICRA`` will require ``python`` v3.6 and ``java`` (we tested in ``openjdk`` 14 2020-03-17).

XICRA requirements
------------------

We include details on the different modules required by ``XICRA`` to successfully run here:

.. toctree:: 
   requirements.rst

.. ########################
.. _install-XICRA:
.. ########################

Installing XICRA
================

If you want to run a stable version of ``XICRA``, the installation can be done with 
``conda`` or ``pip`` packages following instructions in :ref:`Installing from conda<install-from-conda>`
and :ref:`Installing from pip<install-from-pip>`, correspondingly. We encourage you
to install ``XICRA`` and all dependencies using the ``conda`` environment.

..
.. Is it necessary the second part of the installation, From source??
..

On the other hand if you are interested in contributing to ``XICRA`` development, 
running the latest source code, or just like to build everything yourself, you
should follow the :ref:`Installing from source<install-from-source>` instructions. 

Additionally, there are a number of dependencies that might be necessary to 
install or check within your system in both cases. Choose the appopiate choice 
according to your intalling ``XICRA`` option selected.

.. ####################
.. _install-from-conda:
.. ####################

Installing from conda
---------------------

Unfortunately, a couple of executables are not available neither as a ``conda`` or 
``pip`` packages. These packages are ``miraligner`` and ``sRNAbench``. We have generated a
shell script to retrieve ``miraligner`` and include it within your ``conda`` environment.
However, ``sRNAbench`` is no longer available to be downloaded, thus, only the users
with the software already installed will be able to run the ``XICRA`` analysis with 
``sRNAbench``.
To create a new ``conda`` environment, install third party software, install ``XICRA``
and missing dependencies, do as follows:

.. code-block:: sh

   ## clone repo
   git clone https://github.com/HCGB-IGTP/XICRA.git

   ## move to folder
   cd XICRA

   ## create conda environemt
   conda env create -n XICRA -f XICRA_pip/devel/conda/requirements.txt

   ## activate
   conda activate XICRA

   ## install latest python code
   pip install XICRA

   ## install missing software
   sh XICRA_pip/XICRA/config/software/installer.sh
   
To check everything is fine, try executing the ``config`` module:

.. code-block:: sh

   XICRA config
  
We additionally provide a supplementary R package for parsing and 
plotting some ``XICRA`` results, called XICRA.stats_. See additional details here.

Install it in ``R` using:

.. code-block:: sh

   ## Install XICRA.stats version from GitHub:
   ## install.packages("devtools")
   devtools::install_github("HCGB-IGTP/XICRA.stats")
  
.. ##################
.. _install-from-pip:
.. ##################

Installing from pip
-------------------

If you are not using a ``conda`` environment as you might have previously installed all 
dependencies, we encourage you to create a ``python`` environment containing all ``python`` 
modules required for ``XICRA``.  

.. code-block:: sh

   ## create enviroment
   python3 -m venv XICRA_env

   ## activate it
   source XICRA_env/bin/activate

   ## install XICRA and dependencies
   pip install XICRA

   ## execute XICRA
   XICRA -h
  
   
Follow additional `pip installing instructions`_ to learn about installing packages.

Also, as ``XICRA`` relies in multiple dependencies, external packages and third-party 
software, we encourage you to once you install ``XICRA`` using ``pip``. Check for 
dependencies using the ``XICRA config`` module. See
details in ``XICRA config`` module :ref:`section<config-description>`.


.. ##################
.. _install-from-source:
.. ##################

Installing from source
----------------------

Under some circumstancies (develop, bugs fixed, etc) you might be interested 
in obtaining the latest code version. Take into account, that you will need to install 
dependencies and fulfill requirements to have a working distribution. 

.. ##############
.. _get-git-code:
.. ##############

Get source code, for developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have included files in the folder ``devel/conda`` for the build and 
configuration of the ``conda`` package.
The ``XICRA`` project uses git_ as a version control system. To get the code, 
you can grab the latest version from the `XICRA github`_ website, and 
follow the :ref:`Install XICRA from source<install-XICRA-source>` instructions.

Using the command-line, check you have a working distribution of ``git`` by typing 
``git --help`` or install it by typing:

.. code-block:: sh

   sudo apt update
   sudo apt upgrade
   sudo apt install git

Once you have ``git`` installed, to create a new conda environment, install third party
software, install ``XICRA`` and missing dependencies for ``XICRA`` development, do 
as follows:

.. code-block:: sh

   ## clone repo
   git clone https://github.com/HCGB-IGTP/XICRA.git

   ## move to folder XICRA_pip
   cd XICRA/XICRA_pip

   ## create conda environemt
   conda env create -n XICRA -f ./devel/conda/requirements.txt

   ## activate
   conda activate XICRA

   ## install latest python code
   pip install -r ./devel/pypi/requirements.txt
   pip install -e .

   ## install missing software
   sh ./XICRA/config/software/installer.sh
   

.. ###############################
.. _install-XICRA-source:
.. ###############################

Install XICRA from source
^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have ``XICRA`` source code available you only need to include
the ``XICRA`` folder and main script within your path. 

.. code-block:: sh

   export PYTHONPATH=$PYTHONPATH":"$PWD"/XICRA"

   export PATH=$PATH":"$PWD"/XICRA.py"

Take into account that before running ``XICRA`` you have to make sure you have all the 
dependencies fulfilled from section :ref:`Requirements and dependencies<Requirements-dependencies>`.
You can either install them yourself, use appropiate scripts for this purpose or use the ``XICRA config``
module to check, update and install all dependencies required.

.. ###########
.. include:: ../../links.inc