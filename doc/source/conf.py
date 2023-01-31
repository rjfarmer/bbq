# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../src"))


# -- Project information -----------------------------------------------------

project = "BBQ"
copyright = "2023, Rob Farmer"
author = "Rob Farmer"

# The full version, including alpha/beta/rc tags
release = "r23.01.01"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.doctest", "sphinx.ext.napoleon"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

autodoc_mock_imports = []

# standard substitutions
rst_prolog = r"""
.. |MESA I| replace:: `MESA I <https://ui.adsabs.harvard.edu/abs/2011ApJS..192....3P/abstract>`__
.. |MESA II| replace:: `MESA II <https://ui.adsabs.harvard.edu/abs/2013ApJS..208....4P/abstract>`__
.. |MESA III| replace:: `MESA III <https://ui.adsabs.harvard.edu/abs/2015ApJS..220...15P/abstract>`__
.. |MESA IV| replace:: `MESA IV <https://ui.adsabs.harvard.edu/abs/2018ApJS..234...34P/abstract>`__
.. |MESA V| replace:: `MESA V <https://ui.adsabs.harvard.edu/abs/2019ApJS..243...10P/abstract>`__
.. |Msun| replace:: :math:`{\rm M}_\odot`
.. |Lsun| replace:: :math:`{\rm L}_\odot`
.. |Rsun| replace:: :math:`{\rm R}_\odot`
.. |Teff| replace:: :math:`T_{\rm eff}`
.. |logRho| replace:: :math:`\log(\rho/\rm g\,cm^{-3})`
.. |logT| replace:: :math:`\log(T/\rm K)`
.. |chi^2| replace:: :math:`\chi^2`
.. |gpercm3| replace:: :math:`\rm g\,cm^{-3}`
"""

# set default highlighting language
highlight_language = 'fortran'
