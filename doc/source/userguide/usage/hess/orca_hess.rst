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

 * :attr:`~opan.hess.OrcaHess.energy` --
   Electronic energy in :math:`\mathrm{E_h}` (sometimes stored
   spuriously as zero by ORCA)

 * :attr:`~opan.hess.OrcaHess.freqs` --
   Vibrational frequencies in wavenumbers
   :math:`\left(\mathrm{cyc\over cm}\right)`; the six (or five,
   for linear systems) identically-zero frequencies are always at
   the start of the list

 * :attr:`~opan.hess.OrcaHess.hess_path` --
   Full path and filename of the imported HESS file

 * :attr:`~opan.hess.OrcaHess.in_str` --
   Complete contents of the imported HESS file

 * :attr:`~opan.hess.OrcaHess.ir_comps` --
   :math:`x`-, :math:`y`-, and :math:`z`-components in
   :math:`\left(\mathrm{km\over mol}\right)^{1\over 2}` of the transition
   dipole for each normal mode, where row :math:`i` of this matrix matches
   column :math:`i` of :attr:`~opan.hess.OrcaHess.modes`

 * :attr:`~opan.hess.OrcaHess.ir_mags` --
   Squared magnitudes of the transition dipoles
   :math:`\left(\mathrm{T}^2\right)` in
   :math:`\mathrm{km\over mol}`, where element :math:`i` of this array
   matches column :math:`i` of :attr:`~opan.hess.OrcaHess.modes`

 * :attr:`~opan.hess.OrcaHess.joblist` --
   Completion status for the various displacements of a numerical Hessian
   computation, where the entry in the :math:`i^\mathrm{\,th}` row and
   :math:`j^\mathrm{\,th}` column indicates whether the displacement
   of the :math:`i^\mathrm{\,th}` atom along the
   :math:`j`-coordinate has been completed (\ :math:`1`\ ) or not
   (\ :math:`0`\ )

 * :attr:`~opan.hess.OrcaHess.modes` --
   Matrix of normal modes as column vectors, where the full modes matrix
   has been rotation- and translation-purified and mass-weighted, and
   each mode has been separately normalized

 * :attr:`~opan.hess.OrcaHess.mwh_eigvals` --
   1-D array of the eigenvalues of the mass-weighted Hessian, where
   element :math:`i` corresponds to column :math:`i` of both
   :attr:`~opan.hess.OrcaHess.modes` and
   :attr:`~opan.hess.OrcaHess.mwh_eigvecs`

 * :attr:`~opan.hess.OrcaHess.mwh_eigvecs` --
   2-D array of the eigenvectors of the mass-weighted Hessian, where
   column :math:`i` corresponds to element :math:`i` of
   :attr:`~opan.hess.OrcaHess.mwh_eigvals` and column :math:`i` of
   :attr:`~opan.hess.OrcaHess.modes`

 * :attr:`~opan.hess.OrcaHess.num_ats` --
   Number of atoms in the geometry

 * :attr:`~opan.hess.OrcaHess.polders` --
   Derivatives of the elements of the polarizability matrix in units of
   :math:`\mathrm{B}^2~\left(=\mathrm{B^3\over B}\right)`, where
   the :math:`i^\mathrm{\,th}` row contains the derivatives taken with
   respect to the :math:`i^\mathrm{\,th}` coordinate of the geometry, and
   where in each row the matrix element derivatives are presented in the
   in the order of :math:`xx`, :math:`yy`, :math:`zz`,
   :math:`xy`, :math:`xz`, :math:`yz`

 * :attr:`~opan.hess.OrcaHess.raman_acts` --
   Vector of Raman activities in units of
   :math:`\mathrm{\mathring{A}^4 \over u}`

 * :attr:`~opan.hess.OrcaHess.raman_depols` --
   Vector of Raman depolarization ratios

 * :attr:`~opan.hess.OrcaHess.temp` --
   "Actual temperature" reported in the HESS file (sometimes stored
   as a spurious zero value instead of as the actual value used)

Note that not all HESS files contain all of the above data; where data
is absent, in general the respective attribute(s) will be stored as |None|.
In particular:

 * :attr:`~opan.hess.OrcaHess.mwh_eigvals` and
   :attr:`~opan.hess.OrcaHess.mwh_eigvecs` are generally absent unless a
   HESS file is used as input for a mode trajectory run (MTR)

 * :attr:`~opan.hess.OrcaHess.polders`,
   :attr:`~opan.hess.OrcaHess.raman_acts`, and
   :attr:`~opan.hess.OrcaHess.raman_depols` will generally only be present
   if a Raman calculation is requested

 * An exception to this general rule occurs in the case that a numerical
   Hessian is requested but the dipole moment calculation
   is disabled (:code:`%elprop Dipole False end`), where the values of
   :attr:`~opan.hess.OrcaHess.dipders`,
   :attr:`~opan.hess.OrcaHess.ir_comps`, and
   :attr:`~opan.hess.OrcaHess.ir_mags` will be set to zero arrays of the
   appropriate size, instead of |None|

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


