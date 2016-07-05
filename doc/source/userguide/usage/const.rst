.. Usage for const

const -- Constants and Enumerations
===================================

The members of this module fall into three general categories:

 * :ref:`Atomic number/symbol interconversion <usage-const-atnumsym>`

 * :ref:`String and numerical constants <usage-const-constants>`

 * :ref:`Enumerations <usage-const-enums>`


.. _usage-const-atnumsym:

Atomic Numbers & Symbols
------------------------

Conversions between atomic numbers and symbols is provided by two |dict|
members.

:data:`~opan.const.atom_num` provides the atomic number of an
atomic symbol passed in **ALL CAPS**:

>>> opan.const.atom_num['CU']
29

:data:`~opan.const.atom_sym` provides the atomic symbol corresponding to
a given atomic number:

>>> opan.const.atom_sym[96]
'CM'

The elements hydrogen (:math:`Z=1`) through lawrencium (:math:`Z=103`)
are currently supported.

.. _usage-const-constants:

String/Numerical Constants
--------------------------

The purpose of most of these classes and their member values is
sufficiently explained in the respective :ref:`API entries
<api-const-constclasses>`.  The classes anticipated to be most useful to
users are the physical constants of :class:`~opan.const.PHYS`:

.. doctest:: consts

    >>> from opan.const import PHYS
    >>> PHYS.ANG_PER_BOHR
    0.52917721067
    >>> PHYS.LIGHT_SPEED                # Atomic units
    137.036

as well as the string representations of engineering units
of :class:`~opan.const.UNITS`:

.. doctest:: consts

    >>> from opan.const import UNITS
    >>> from opan.const import EnumUnitsRotConst as EURC
    >>> UNITS.rot_const[EURC.INV_INERTIA]
    '1/(amu*B^2)'
    >>> UNITS.rot_const[EURC.ANGFREQ_SECS]
    '1/s'

.. testcleanup:: consts

    del PHYS, UNITS, EURC

Two of the remaining classes, :class:`~opan.const.DEF` and
:class:`~opan.const.PRM`, define default values that are primarily relevant
to programmatic use of ``opan``. In unusual circumstances, though, they may
be useful in console interactions.

:class:`~opan.const.CIC` currently covers a very limited scope (the minimum
and maximum atomic numbers implemented) and will likely not be useful
at the console.

.. _usage-const-enums:

Enumerations
------------

From the perspective of the end user, enumerations in ``opan`` are
"functional types," which don't need instantiation before use:

>>> opan.const.EnumDispDirection.NO_DISP
'NO_DISP'

The enumeration values are simple strings:

>>> type(opan.const.EnumDispDirection.NO_DISP)
<class 'str'>

While this implementation is susceptible to accidental mixing of enumerated
types, it has the advantage of allowing simple |str| inputs to functions
expecting enumerated values.  This is anticipated to be useful in
console-level interactions with a variety of program elements.  For example,
the engineering units to be output from :func:`opan.utils.inertia.rot_consts`
can be specified simply with the appropriate string, instead of the
fully specified enumeration object:

.. testsetup:: enums

    import opan, numpy as np
    some_geom = np.array([2, 7, -3, -2, 15, -5])
    some_masses = np.array([5, 12])

.. doctest:: enums

    >>> from opan.utils.inertia import rot_consts
    >>> from opan.const import EnumUnitsRotConst as EURC
    >>> rot_consts(some_geom, some_masses, 'INV_INERTIA')       # Works fine
    array(...)
    >>> rot_consts(some_geom, some_masses, EURC.INV_INERTIA)     # Also works
    array(...)

As noted in the API documentation for :class:`~opan.const.EnumIterMeta`,
both iteration and membership testing with
":ref:`in <python:in>`" are supported:

.. doctest:: enums

    >>> 'NO_DISP' in opan.const.EnumDispDirection
    True
    >>> [e for e in sorted(opan.const.EnumDispDirection)]
    ['NEGATIVE', 'NO_DISP', 'POSITIVE']

.. testcleanup:: enums

    del some_geom, some_masses







