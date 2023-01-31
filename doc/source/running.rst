Running
=======

When ``./bbq`` is invoked it will search for a file called ``inlist`` in the current working directory. 
You can override this by providing a path to the inlist on the command line ``./bbq /path/to/other/inlist``.

An empty ``inlist`` file is provide in the root folder. See the ``defaults`` folder for the available options and their default options.


Modes
~~~~~


There are currently three modes that can be ran, each mode has its own inlist options:

* Sampling
* Random Sampling
* Profile

A note on inputs/outputs, the base units for all quantities are: temperatures in Kelvin, densities in g/cm^{-3}, time in seconds, and composition in abundances. 
Some quantities may be in log10 when specified.

The current set of default options can be found in the ``defaults/`` folder. Like MESA if you change anythin in the ``defaults`` folder you must recompile
``bbq`` for the changes to take effect.


Sampling
--------

In this mode we read an input file line by line and process each line as an independent set of inputs. Inputs are specified in the ``sampling`` inlist

Random
------

Similiar to ``sampling`` except the inputs are drawn from log flat distributions and processed as an independent set of inputs. Inputs are specified in the ``random`` inlist


Profile
-------

Here we read one set of initial compositions and then follow that material as it moves along the thermodyanmic trajectory given by the input file.
Inputs are given by the ``profile`` inlist.

Hydrostatic
-----------

Similar to profile except we maintain a constant T/Rho during the burn.
