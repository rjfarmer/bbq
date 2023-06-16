.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7585202.svg
   :target: https://doi.org/10.5281/zenodo.7585202
.. image:: https://readthedocs.org/projects/stellar-bbq/badge/?version=latest
    :target: https://stellar-bbq.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


BBQ
===

``bbq`` is a one-zone nuclear network solver that uses MESA's microphysics for solving the 
change in composition from nuclear burning. It does this by integrating the change in 
composition due to nuclear burning at constant temperature and density with a semi-implicit midpoint rule.
More details about the MESA implementation can be found in `Section 10.2 of Jermyn et al 2023 <https://ui.adsabs.harvard.edu/abs/2022arXiv220803651J/abstract>`_.


Installation
------------

First ensure that you have installed `MESA <https://docs.mesastar.org/en/release-r22.11.1/installation.html>`_ and have the
environment variable ``MESA_DIR`` pointing to your MESA installation.

Currently supported MESA versions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* 22.11.01
* 23.05.1

Then run::

    ./mk

You will now have a executable called ``bbq`` which can be ran with ``./bbq`` to run the code. 
This file can be moved outside of the 
source folder if wanted, it has no dependencies on the src folder.

Note that the test suite will likely only pass with the newest MESA version supported. This is due to changes in the reaction
rates between MESA versions which mean the test suite output needs updating for each release.


Documentation
~~~~~~~~~~~~~

Online documentation can be found `readthdocs <https://stellar-bbq.readthedocs.io/en/latest/>`_.

To generate the documentation locally::

    cd doc
    pip install -r requirments.txt
    make html


Running
-------

When ``./bbq`` is invoked it will search for a file called ``inlist`` in the current working directory. 
You can override this by providing a path to the inlist on the command line ``./bbq /path/to/other/inlist``.


Modes of operation
------------------

There are currently four modes that can be ran, each mode has its own inlist options:

* Sampling
* Random Sampling
* Profile
* hydrostatic


A note on inputs/outputs, the base units for all quantities are: temperatures in Kelvin, densities in g/cm^{-3}, time in seconds, and composition in mass fractions. 
Some quantities may be in log10 when specified.

The current set of default options can be found in the ``defaults/`` folder. Like MESA if you change anything in the ``defaults`` folder you must recompile
``bbq`` for the changes to take effect.


Sampling
~~~~~~~~

In this mode we read an input file line by line and process each line as an independent set of inputs. Inputs are specified in the ``sampling`` inlist

Random
~~~~~~

Similiar to ``sampling`` except the inputs are drawn from log flat distributions and processed as an independent set of inputs. Inputs are specified in the ``random`` inlist


Profile
~~~~~~~

Here we read one set of initial compositions and then follow that material as it moves along the thermodyanmic trajectory given by the input file.
Inputs are given by the ``profile`` inlist.

Hydrostatic
~~~~~~~~~~~

Similar to profile except we maintain a constant T/Rho during the burn. Inputs are given by the ``hydrostatic`` inlist.


Other inlists
-------------

Other inlists exist which specify the nuclear physics and EOS used.

BBQ inlist
~~~~~~~~~~

The ``bbq`` inlist is the main driver inlist and is where you set which mode to operate in as well as the nuclear network used.

eos inlist
~~~~~~~~~~

The ``eos`` inlist is the same as the MESA EOS inlist and can contain anything that inlist specifies. It is the only inlist that allows
for nesting of other inlists.
See `MESA's eos options <https://docs.mesastar.org/en/release-r22.11.1/reference/eos.html>`_ for the full set of supported options.


nuclear inlist
~~~~~~~~~~~~~~

The ``nuclear`` inliust contains many of the nuclear related inlist options from MESA's `star_job <https://docs.mesastar.org/en/release-r22.11.1//reference/star_job.html>`_ and `controls <https://docs.mesastar.org/en/release-r22.11.1/reference/controls.html>`_ inlists. 


Citing this work
----------------

See the CITATIONS.bib file for a set of bibtex entries to cite. At a minimum you should cite this repository via via `Zenodo <https://doi.org/10.5281/zenodo.7585201>`_ with the version used, as well as 
`Jermyn et al 2023 <https://ui.adsabs.harvard.edu/abs/2022arXiv220803651J/abstract>`_. You should also cite the source of the microphysics you use (reaction rates, screening, equation of state, etc).
The `MESA CITATIONS.bib <https://github.com/MESAHub/mesa/blob/main/CITATIONS.bib>`_ file contains many of these references, but it is up to you to ensure they are correct. 
