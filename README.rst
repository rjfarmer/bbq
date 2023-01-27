bbq
===

`bbq` is a one-zone nuclear network solver that uses MESA's microphysics for solving the change in composition from nuclear burning.
It does this by integrating the change in composition due to nuclear burning at constant temperature and density with a semi-implicit midpoint rule.
More details about the MESA implementation can be found in `Section 10.2 <https://ui.adsabs.harvard.edu/abs/2022arXiv220803651J/abstract>`_.


Installation
------------

First ensure that you have installed `MESA <https://docs.mesastar.org/en/release-r22.05.1/installation.html>`_ and have the
environment variable `MESA_DIR` pointing to your MESA installation.

Currently supported MESA versions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* 22.05.01
* 22.11.01

Previous versions may work but depend on how the MESA internals have changed.


Then run::

    ./mk

You will now have a executable called `bbq` which can be ran with `./bbq` to run the code. This file can be moved outside of the 
source folder if wanted, it has no dependencies on the src folder.


Running
-------

When `./bbq` is invoked it will search for a file called `inlist` in the current working directory. You can over this by providing path to the inlist on the command line `./bbq /path/to/other/inlist`.


Modes of operation
------------------

There are currently three modes that can be ran, each mode has its own inlist options:

* Sampling
* Random Sampling
* Profile

In the following section we describe each mode and its inlist options

A note on inputs/outputs, the base units for all quantities are: temperatures in Kelvin, densities in g/cm^{-3}, time in seconds, and composition in abundances. Some quantities may be in log base 10 when specified.


bbq
~~~

Common options for all modes::

    &bbq

      ! Nuclear network in MESA format. If path is not fully specified then we look in $MESA_DIR/data/net_data/nets
      net_name = 'approx21.net'

      ! Screening mode
      screening_mode = 'chugunov'

      ! Scale factor for weak reactions
      weak_rate_factor = 1

      ! Solver tolerances
      max_steps = 1000000
      eps = 1d-8
      odescal = 1d-10
      stptry = 0


      ! What mode to run, must pick one
      use_input_file = .false. ! Simply read from file and output burn data
      use_random_sampling =.false. ! Randomly sample inputs
      use_profile=.false. ! Read a profile of (T,Rho) and burn material

      ! Only outputs the list of isotopes then stop
      just_write_isos=.false.

      ! whether to output the list of isotopes in net_iso order
      write_iso_list = .false.
      iso_list_filename = 'isos.txt'

      ! How often to force the output to disk
      flush_freq = 50 

    /



For specifying the composition in input files (and the order in the output files) this must always be in MESA's `net_iso` order. To determine the `net_iso` order for a nuclear network set `just_write_isos=.true.` and `iso_list_filename` and run `bbq`. This will output in the `iso_list_filename` file the isotopes in `net_iso` order, one per line.



Sampling
~~~~~~~~

In this mode we read the data from an input file and process each line one by one, assuming that each line is independent of each other::

    &sampling

        ! File containing input data as: 
        ! log(time/s) log(T/K) log(rho/gcm^{-3}) composition
        input_filename = ''

        ! Outputs:
        !eps_nuc(erg/g/gs) eps_neu(erg/g/gs) composition
        output_filename = ''

    /

Random Sampling
~~~~~~~~~~~~~~~

In this mode we randomly sample the input options::

    &random

      ! Where to output starting points for the sample
      ! log(time/s) log(T/K) log(rho/gcm^{-3}) composition
      output_starting_filename 

      ! Where to output final composition of step
      !eps_nuc(erg/g/gs) eps_neu(erg/g/gs) composition
      output_ending_filename 

      ! How many samples to draw, negative means unlimited
      num_samples 

      ! Min and Max values for the log(time/s) to integrate for
      log_time_min  
      log_time_max 

      ! Min and Max values for the log(T/Ks) temperature
      log_temp_min 
      log_temp_max 

      ! Min and Max values for the log(rho/gcm^{-3}) density
      log_rho_min 
      log_rho_max 

      ! Min and Max values for the log(Xa) abundance
      log_xa_min 
      log_xa_max 

      ! Place a limit on the abundance of free neutrons and protons
      neut_prot_limit_frac 

      ! Random seed, if negative use a different one each time 
      seed = 42 

    /


Profile
~~~~~~~

In this mode we read the data from an input file and process each line but assume that the composition is moving between each set of (time,Temp,rho) triplets::

    &profile
        ! Input file with
        ! log(time/s) log(T/K) log(rho/gcm^{-3})
        input_filename 

        ! Input file the abundances in net_iso order (one iso per line)
        input_composition_filename

        ! Outputs
        ! total_time(s) dt(s) log(T/K) log(rho/gcm^{-3})
        output_filename 

        ! After processing input_filename what to do at the end,
        ! if this is false we just stop
        ! if this is true we will repeat the thermodynamic trajectory but in reverse order and do this num_loops number of times.
        ! This can be thought of if a packet of material was being burnt while being covectively mixed and we wish to follow the flow as it rises and falls multiple times.
        reflective_boundaries=.true.
        num_loops = 1

    / 



eos
~~~

This isn't a standalone mode but just specifies the EOS choices (this is needed for the elctron degeneracy parameter that the weak rates need)
This is exactly the same as MESA's normal eos inlist and supports nested calls to other eos inlists.

See `MESA's eos options <https://docs.mesastar.org/en/release-r22.05.1/reference/eos.html>`_ for the full set of supported options.
