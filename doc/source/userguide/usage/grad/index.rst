.. Usage for grad core/base/superclass

grad
====

While not as important as Hessians for VPT2 calculations,
nuclear gradients can be used to estimate the significance of certain
terms in the VPT2 expansion prior to calculating the associated Hessians,
potentially reducing the associated computational demand [Bar11]_. Other
calculations envisioned for eventual Open Anharmonic
implementation would also make more direct use of gradient data.

In order to enable automated parsing of gradient data from various
computational software packages, a common interface must be defined
by which the automated components of ``opan`` can access needed data in
a uniform manner.  This interface is defined by the abstract
:class:`~opan.grad.SuperOpanGrad` superclass, which requires that the
following members be defined for all of its subclasses:

 * :attr:`gradient` -- 1-D |nparray| of |npfloat| -- Gradient data in
   :math:`\left(\mathrm{E_h\over B}\right)` |units|

 * :attr:`geom` -- 1-D |nparray| of |npfloat| -- Geometry data in
   :math:`\mathrm B` |units|

 * :attr:`atom_syms` -- |list| of |str| -- Atomic symbols in **ALL CAPS**

Beyond these required members, the below subpages describe
the various imported data that are available in the subclasses of
:class:`~opan.grad.SuperOpanGrad`.

**Implemented Subclasses**

.. toctree::
    :maxdepth: 1

    orca_engrad

