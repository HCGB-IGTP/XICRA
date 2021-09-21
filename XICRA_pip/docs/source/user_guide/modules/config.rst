.. ############################
.. _config-description:

config
======

The ``XICRA`` pipeline provides this module ``config`` to help the user
to configurate the multiple dependencies and requirements.

We encourage ``XICRA`` users to run this module after the installation to
check whether the multiple requirements and dependencies are fulfilled.

.. seealso:: Additional information on ``XICRA`` configuration and requirements

   - :doc:`Installation <../installation/installing>` 
   
   - :doc:`Requirements <../installation/requirements>` 
   
   - :doc:`config module (API) <../../api/modules/config>`


How to run the config module
----------------------------

Once you have installed ``XICRA`` you should be able to run in the command line the pipeline.

If you type ``XICRA`` you should see a prompt with the different modules available. Following 
the pipeline name type the module of interest, in this case, config. As an example:

.. code-block:: sh

   XICRA config -h

The different options and parameters for this module should appear in the command line prompt. Here we summarized them:


.. function:: Module XICRA config
   
   :param -h, --help: Show this help message and exit. 
   :param --debug: Show additional message for debugging purposes.

Basically, this ``XICRA config`` module allows the user to check if the requirements are fulfilled. 

.. code-block:: sh

   XICRA config

.. include:: ../../links.inc