Installation
============

First ensure that you have installed `MESA <https://docs.mesastar.org/en/release-r22.11.1/installation.html>`_ and have the
environment variable ``MESA_DIR`` pointing to your MESA installation.

Currently supported MESA versions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* 22.11.01


Then run::

    git clone git@github.com:rjfarmer/bbq.git
    cd bbq
    ./mk

You will now have a executable called ``bbq`` which can be ran with ``./bbq`` to run the code. 
This file can be moved outside of the 
source folder if wanted, it has no dependencies on the src folder.
