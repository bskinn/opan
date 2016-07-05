.. Usage for hess core/base/superclass

hess -- Hessian Data
====================

Nuclear Hessians are key for VPT2 calculations, as they provide the third-
and fourth-order force constants for a system via finite difference
differentiation of the second-order constants calculated by standard
harmonic-oscillator methods [Bar11]_ [Blo12]_.

.. todo::

    Link to HO writeup once composed

In order to enable automated parsing of Hessian data from various
computational software packages, a common interface must be defined
by which the automated components of ``opan`` can access needed data in
a uniform manner.  This interface is defined by the abstract
:class:`~opan.hess.SuperOpanHess` superclass, which requires that the
following members be defined for all of its subclasses:

 * :attr:`hess` -- 2-D |nparray| of |npfloat| -- Hessian data in
   :math:`\left(\mathrm{E_h\over B^2}\right)` |units|

 * :attr:`geom` -- 1-D |nparray| of |npfloat| -- Geometry data in
   :math:`\mathrm B` |units|

 * :attr:`atom_syms` -- |list| of |str| -- Atomic symbols in **ALL CAPS**

Beyond these required members, the below subpages describe
the various imported data that are available in the subclasses of
:class:`~opan.hess.SuperOpanHess`.

**Implemented Subclasses**

.. toctree::
    :maxdepth: 1

    orca_hess

