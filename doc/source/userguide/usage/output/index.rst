.. Usage for output

output -- Computation Outputs
=============================

The top-level structure required for computation outputs is not yet fully
defined.  Once the particulars of the VPT2 computations are implemented,
a ``SuperOpanOutput`` superclass will likely be created to enforce
this structure, as with :class:`~opan.grad.SuperOpanGrad` and
:class:`~opan.hess.SuperOpanHess`.

At this time, preliminary support for parsing |orca| output is provided
via :class:`~opan.output.OrcaOutput`.

.. toctree::
    :maxdepth: 1

    orca_output


