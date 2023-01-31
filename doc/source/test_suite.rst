Test suite
==========

The folder ``test`` contains a number of example test cases.

All tests can be run by invoking ``./each_test_run`` which will cycle over each test and compare its output with a reference file. The script ``./each_test_update`` will update the reference file for each test case.


Anatomy of a test case
~~~~~~~~~~~~~~~~~~~~~~

A new test case can be built by first copying the folder ``sample_test_suite`` to a new folder. Inside this folder are a number of scripts ``rn``, ``test``, ``update``, ``ck`` and a basic inlist ``inlist``.

The ``rn`` and ``test`` scripts should not normally need to be modified. However the ``update`` and ``ck`` scripts should be updated when making a new test case to ensure that the appropriate files are being tested. Generally we want to run ``diff`` comparing the output to a reference file, but sometimes there may be multiple files that need testing.


Test cases
~~~~~~~~~~

big_net
-------

This tests loading a running a large 495 isotope nuclear network. For the first time after a clean ``MESA`` installation is may take some time to rebuild the rate cache.


isos
----

A short test just to write out the network in ``net_iso`` order.

ni56
----

Runs a custom nuclear network that follows Ni56 -> Co56 -> Fe56. While it runs in sample mode (so each input line is treated independently) the dt is increasing so it acts like you are following a thermodyanmic trajectory.

profile
-------

Follows the composition over a thermodyanmic trajectory, where the composition is set once at the start of the run and evolved forward in time given the (logT, logRho).

random
------

Generate a run in sample mode, where the inputs are selected randomly. Note the use of a random seed, if that is set to >0 to ensure reproducible results, if negtaive the see is set based on the clock.

sample
------

Run in sample mode (so each input line is treated independently).
