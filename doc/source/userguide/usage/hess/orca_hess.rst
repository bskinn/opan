.. Usage for OrcaHess

.. testsetup:: import

    import opan, os
    olddir = os.getcwd()
    os.chdir(os.path.join('source', '_static'))

.. testsetup:: usage

    import opan, os
    olddir = os.getcwd()
    os.chdir(os.path.join('source', '_static'))
    h = opan.hess.OrcaHess(path='h2o.hess')

.. testcleanup:: *

    os.chdir(olddir)

OrcaHess
========

In addition to :attr:`~opan.hess.OrcaHess.hess`,
:attr:`~opan.hess.OrcaHess.geom`, and
:attr:`~opan.hess.OrcaHess.atom_syms`
as required by the :class:`~opan.hess.SuperOpanHess`
:doc:`specification <index>`, |orca| HESS files expose the following
attributes:

 * :attr:`~opan.hess.OrcaHess.atom_masses` --
   Atom masses as reported in the geometry block, in :math:`\mathrm u`

 * :attr:`~opan.hess.OrcaHess.dipders` --
   Dipole derivatives in :math:`\mathrm{D_a\over B}` atomic |units|,
   where the element in the :math:`i^\mathrm{\,th}` row and
   :math:`j^\mathrm{\,th}` column is
   the derivative of the :math:`j`-component of the dipole moment with
   respect to the :math:`i^\mathrm{\,th}` coordinate of the geometry

 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --


 * :attr:`~opan.hess.OrcaHess.` --

Note that not all HESS files contain all of the above data; where data
is absent, in general the respective attribute(s) will be stored as |None|.
For certain data which are expected to reside in *all* HESS files
(those annotated as *(required)* in the
:ref:`instance variables list <hess-orcahess-instancevars>` for
:class:`~opan.hess.OrcaHess`), a :class:`~opan.error.HessError` will be
raised if any are absent.

The public class :class:`OrcaHess.Pat <opan.hess.OrcaHess.Pat>` contains
|re.compile| patterns used during file import. Their usefulness thus may be
limited.

Import a HESS file by passing its full path and name to the
:class:`~opan.hess.OrcaHess` constructor via the `path` keyword argument:

.. doctest:: import

    h = opan.hess.OrcaHess(path='h2o.hess')

The contents of the file are accessible as simple attributes:

.. doctest:: usage

    >>> h.hess[0:6,0:6]
    array([[  5.48742000e-01,  -1.01850000e-01,  -1.00000000e-06,
             -5.00222000e-01,   2.25860000e-02,   1.00000000e-06],
           [ -1.01817000e-01,   5.48777000e-01,  -2.00000000e-06,
              7.92730000e-02,  -4.85500000e-02,  -0.00000000e+00],
           [ -1.28000000e-04,   1.03000000e-04,   1.35300000e-03,
              1.49000000e-04,   2.50000000e-05,  -6.76000000e-04],
           [ -5.00624000e-01,   7.90960000e-02,   0.00000000e+00,
              5.09553000e-01,  -7.24600000e-02,  -0.00000000e+00],
           [  2.27460000e-02,  -4.84520000e-02,   1.00000000e-06,
             -7.27740000e-02,   5.77420000e-02,  -1.00000000e-06],
           [ -9.50000000e-05,  -1.99000000e-04,  -7.27000000e-04,
              9.00000000e-06,   6.10000000e-05,   5.10000000e-04]])
    >>> h.geom
    array([-0.088833, -0.088832,  0.      ,  1.721059, -0.31111 ,  0.      ,
           -0.311106,  1.721053,  0.      ])
    >>> h.atom_syms
    ['O', 'H', 'H']
    >>> h.num_ats
    3
    >>> h.freqs
    array([    0.      ,     0.      ,     0.      ,     0.      ,
               0.      ,     0.      ,  1610.279974,  3761.722714,  3848.311829])

Again, if data is not available it will be stored as |None|:

.. doctest:: usage

    >>> h.polders is None
    True


