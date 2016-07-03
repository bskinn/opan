.. Usage for OrcaEngrad

.. testsetup:: import

    import opan, os
    olddir = os.getcwd()
    os.chdir(os.path.join('source', '_static'))

.. testsetup:: usage

    import opan, os
    olddir = os.getcwd()
    os.chdir(os.path.join('source', '_static'))
    oe = opan.grad.OrcaEngrad(path='h2o.engrad')

.. testcleanup:: *

    os.chdir(olddir)

OrcaEngrad
==========

In addition to :attr:`~opan.grad.OrcaEngrad.gradient`,
:attr:`~opan.grad.OrcaEngrad.geom`, and
:attr:`~opan.grad.OrcaEngrad.atom_syms`
as required by the :class:`~opan.grad.SuperOpanGrad`
:doc:`specification <index>`, |orca| ENGRAD files expose the following
attributes:

 * :attr:`~opan.grad.OrcaEngrad.energy` --
   Single-point energy for the related geometry, in
   :math:`\mathrm{E_h}`

 * :attr:`~opan.grad.OrcaEngrad.in_str` --
   Full contents of the ENGRAD file

 * :attr:`~opan.grad.OrcaEngrad.num_ats` --
   Number of atoms in the system

The public class :class:`OrcaEngrad.Pat <opan.grad.OrcaEngrad.Pat>` contains
|re.compile| patterns used during file import. Their usefulness thus may be
limited.

Import an ENGRAD file by passing its full path and name to the
:class:`~opan.grad.OrcaEngrad` constructor via the `path` keyword argument:

.. doctest:: import

    oe = opan.grad.OrcaEngrad(path='h2o.engrad')

The contents of the file are accessible as simple attributes:

.. doctest:: usage

    >>> oe.gradient
    array([ -2.33839000e-07,  -2.33870000e-07,  -8.00000000e-12,
             2.26705000e-07,   7.13300000e-09,   2.00000000e-12,
             7.13400000e-09,   2.26737000e-07,   5.00000000e-12])
    >>> oe.geom
    array([ 0.1613726,  0.1613726,  0.       ,  1.9824472, -0.0651211,
            0.       , -0.0651211,  1.9824472,  0.       ])
    >>> oe.atom_syms
    ['O', 'H', 'H']
    >>> oe.num_ats
    3

This ENGRAD file is for an optimized geometry of water, and thus all elements
of the gradient are quite small.


