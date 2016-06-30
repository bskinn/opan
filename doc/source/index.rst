.. opan documentation master file, created by
   sphinx-quickstart on Sun Sep 20 15:06:37 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to Open Anharmonic!
===========================

Open Anharmonic is a Python wrapper for computational chemistry
software packages intended to enable VPT2 computation of anharmonic
vibrational constants. The code is still in the preliminary stages of
development; no VPT2 functionality is yet available.

Other types of calculations are under consideration.

An adjunct goal of the project is to expose an API providing
convenient access to various results of standalone calculations, as well
as tools to manipulate those results.
In particular, :class:`~opan.xyz.OpanXYZ` and the subclasses of
:class:`~opan.grad.SuperOpanGrad` and
:class:`~opan.hess.SuperOpanHess` are anticipated to be particularly
useful.

Due to the large number of possible variations of computational
runs, parsing of output files is challenging, and only a small number of
run types have been `implemented to date <api/output.html>`__. More are planned,
but are currently low priority.

Open Anharmonic is available on PyPI as ``opan``::

    pip install opan

See the :doc:`dependencies page <dependencies>` for package dependencies
and compatible versions. Note that due to common complications on Windows
systems, dependencies are **NOT** set to automatically install.

The source repository for Open Anharmonic can be found at:

    https://www.github.com/bskinn/opan

Bug reports and feature requests can be submitted as GitHub Issues. Other
feedback is welcomed at::

    bskinn at alum dot mit dot edu

**Contents**

.. toctree::
    :maxdepth: 1

    User's Guide <userguide/index>
    API <api/index>
    conventions
    software
    Dependencies <dependencies>
    units
    references
    doc_todo


**Indices and Tables**

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


