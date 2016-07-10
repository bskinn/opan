#-------------------------------------------------------------------------------
# Name:        hess
# Purpose:     Encapsulation of parsing/handling of data from Hessian files.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     14 Nov 2014
# Copyright:   (c) Brian Skinn 2016
# License:     The MIT License; see "license.txt" for full license terms
#                   and contributor agreement.
#
#       This file is part of opan (Open Anharmonic), a system for automated
#       computation of anharmonic properties of molecular systems via wrapper
#       calls to computational/quantum chemical software packages.
#
#       http://www.github.com/bskinn/opan
#
#-------------------------------------------------------------------------------

""" Module implementing imports of Hessian data from external computations.

The abstract superclass :class:`SuperOpanHess` defines a common initializer
and common method(s)  for use by subclasses importing
Hessian data from external computational packages.

|

**Implemented Subclasses**

:class:`OrcaHess` -- Imports '.hess' files from |orca|

|

**Requirements**

 *  The import for each external software package SHOULD have
    its own subclass.

 *  Each subclass MUST implement a ``_load(**kwargs)`` method as the
    entry point for import of Hessian data.

 *  The Hessian data MUST be stored:

    *   In the instance member ``self.hess``

    *   As a two-dimensional |nparray|

    *   With `dtype` descended from |npfloat|

    *   In |units| of Hartrees per Bohr-squared
        :math:`\\left(\\frac{\\mathrm{E_h}}{\\mathrm B^2}\\right)`

    *   With elements arranged as:

    .. math::

        \\left[
        \\begin{array}{ccccccc}
            \\frac{\\partial^2 E}{\\partial x_1^2} &
                \\frac{\\partial^2 E}{\\partial x_1\\partial y_1} &
                \\frac{\\partial^2 E}{\\partial x_1\\partial z_1} &
                \\frac{\\partial^2 E}{\\partial x_1\\partial x_2} &
                \\dots &
                \\frac{\\partial^2 E}{\\partial x_1\\partial y_N} &
                \\frac{\\partial^2 E}{\\partial x_1\\partial z_N} \\\\
            \\frac{\\partial^2 E}{\\partial y_1\\partial x_1} &
                \\frac{\\partial^2 E}{\\partial y_1^2} &
                \\frac{\\partial^2 E}{\\partial y_1\\partial z_1} &
                \\frac{\\partial^2 E}{\\partial y_1\\partial x_2} &
                \\dots &
                \\frac{\\partial^2 E}{\\partial y_1\\partial y_N} &
                \\frac{\\partial^2 E}{\\partial y_1\\partial z_N} \\\\
            \\frac{\\partial^2 E}{\\partial z_1\\partial x_1} &
                \\frac{\\partial^2 E}{\\partial z_1\\partial y_1} &
                \\frac{\\partial^2 E}{\\partial z_1^2} &
                \\frac{\\partial^2 E}{\\partial z_1\\partial x_2} &
                \\dots &
                \\frac{\\partial^2 E}{\\partial z_1\\partial y_N} &
                \\frac{\\partial^2 E}{\\partial z_1\\partial z_N} \\\\
            \\frac{\\partial^2 E}{\\partial x_2\\partial x_1} &
                \\frac{\\partial^2 E}{\\partial x_2\\partial y_1} &
                \\frac{\\partial^2 E}{\\partial x_2\\partial z_1} &
                \\frac{\\partial^2 E}{\\partial x_2^2} &
                \\dots &
                \\frac{\\partial^2 E}{\\partial x_2\\partial y_N} &
                \\frac{\\partial^2 E}{\\partial x_2\\partial z_N} \\\\
            \\vdots & \\vdots & \\vdots & \\vdots & \\ddots &
                \\vdots & \\vdots \\\\
            \\frac{\\partial^2 E}{\\partial y_N\\partial x_1} &
                \\frac{\\partial^2 E}{\\partial y_N\\partial y_1} &
                \\frac{\\partial^2 E}{\\partial y_N\\partial z_1} &
                \\frac{\\partial^2 E}{\\partial y_N\\partial x_2} &
                \\dots &
                \\frac{\\partial^2 E}{\\partial y_N^2} &
                \\frac{\\partial^2 E}{\\partial y_N\\partial z_N} \\\\
            \\frac{\\partial^2 E}{\\partial z_N\\partial x_1} &
                \\frac{\\partial^2 E}{\\partial z_N\\partial y_1} &
                \\frac{\\partial^2 E}{\\partial z_N\\partial z_1} &
                \\frac{\\partial^2 E}{\\partial z_N\\partial x_2} &
                \\dots &
                \\frac{\\partial^2 E}{\\partial z_N\\partial y_N} &
                \\frac{\\partial^2 E}{\\partial z_N^2} \\\\
        \\end{array}
        \\right]

    .. note::

        The Hessian is elsewhere assumed to be symmetric
        (real-Hermitian), and thus MUST be returned as such here.
        Symmetric character is **NOT** explicitly checked, however!

 *  The geometry data MUST be stored:

    *   In the instance member ``self.geom``

    *   As a one-dimensional |nparray|

    *   With `dtype` descended from |npfloat|

    *   In |units| of Bohrs :math:`\\left(\\mathrm B\\right)`

    *   With elements ordered as:

    .. math::

        \\left[
        \\begin{array}{cccccccc}
            x_1 & y_1 & z_1 & x_2 & y_2 & \\dots & y_N & z_N \\\\
        \\end{array}
        \\right]

 *  The atoms list MUST be stored:

    *   In the instance member ``self.atom_syms``

    *   As a |list| of |str|, with each atom specified by an **all-caps**
        atomic symbol (:data:`opan.const.atom_sym` may be helpful)

 *  Subclasses MAY define an unlimited number of methods,
    class variables, and/or
    instance variables in addition to those defined above, of
    unrestricted type.

|

**Superclass**

.. autoclass:: SuperOpanHess

|

**Subclasses**

.. autoclass:: OrcaHess(path='...')


"""


# Imports


# Debug constant
_DEBUG = False



class SuperOpanHess(object):
    """ Abstract superclass of Hessian import classes.

    Performs the following actions:

    1.  Ensures that the abstract superclass is not being instantiated,
        but instead a subclass.
    2.  Calls the ``_load()`` method on the subclass, passing all `kwargs`
        through unmodified.
    3.  Typechecks the ``self.hess``, ``self.geom``, and
        ``self.atom_syms`` required members for existence, proper data type,
        and properly matched lengths/dimensions.

    The checks performed in step 3 are primarily for design-time member
    and type enforcement during development of subclasses for new
    external software packages, rather than run-time data checking.
    It is RECOMMENDED to include robust data validity checking inside
    each subclass, rather than relying on these tests.

    |

    **Methods**

    .. automethod:: check_geom(coords, atoms[, tol])

    """

    # Imports
    from .const import DEF as _DEF

    def __init__(self, **kwargs):
        """ Conformity-checking initalizer for Hessian data loader classes.
        """

        # Imports
        from .error import HessError as HErr
        from .utils import assert_npfloatarray as a_npfa
        from .const import atom_num
        import numpy as np

        # Check for abstract base class
        if type(self) == SuperOpanHess:
            raise NotImplementedError("SuperOpanHess base class is abstract")
        ## end if

        # Call the subclass _load method, passing in the keyword arguments
        #  wholesale
        self._load(**kwargs)

        # Define common error source string
        srcstr = "{0} with args {1}".format(self.__class__, str(kwargs))

        # All proofreading here is excluded from coverage because
        #  all of these should be handled at the subclass level.

        # Proofread types for Hessian and geometry
        a_npfa(self, 'hess', 'Hessian', HErr, HErr.BADHESS, srcstr)
        a_npfa(self, 'geom', 'geometry', HErr, HErr.BADGEOM, srcstr)

        # Hessian and geometry shape
        if (len(self.hess.shape) != 2 or
                            self.hess.shape[0] % 3 != 0 or
                            self.hess.shape[0] != self.hess.shape[1]):
            raise HErr(HErr.BADHESS, # pragma: no cover
                        "Hessian is not a 3N x 3N matrix", srcstr)
        ## end if
        if self.geom.shape != (self.hess.shape[0],):
            raise HErr(HErr.BADGEOM, # pragma: no cover
                        "Geometry shape does not match Hessian shape", srcstr)

        # Ensure atomic symbols length matches & they're all valid
        if not hasattr(self, 'atom_syms'): # pragma: no cover
            raise HErr(HErr.BADATOM, "Atoms list not found", srcstr)
        ## end if
        if type(self.atom_syms) is not type([1,2,3]): # pragma: no cover
            raise HErr(HErr.BADATOM, "Atoms list is not a list", srcstr)
        ## end if
        if 3*len(self.atom_syms) != self.hess.shape[0]: # pragma: no cover
            raise HErr(HErr.BADATOM, "Atoms list is not length-N", srcstr)
        ## end if
        if not all(map(lambda v: v in atom_num, self.atom_syms)):
            raise HErr(HErr.BADATOM,    # pragma: no cover
                    "Invalid atoms in list: {0}".format(self.atom_syms),
                    srcstr)
        ## end if

    ## end def __init__


    def check_geom(self, coords, atoms, tol=_DEF.HESS_COORD_MATCH_TOL):
        """ Check for consistency of Hessian geometry with input coords/atoms.

        The cartesian coordinates associated with a Hessian object
        are considered consistent with the input `coords` and `atoms`
        if each component matches to within `tol` and all atoms
        are identical.  If
        `coords` or `atoms` vectors are passed that are of different
        length than those stored in the instance, a |False| value is
        returned, rather than an exception raised.

        Parameters
        ----------
        coords
            length-3N |npfloat_| --
            Vector of stacked 'lab-frame' Cartesian coordinates

        atoms
            length-N |str| or |int| --
            Vector of atom symbols or atomic numbers

        tol
            |float|, optional --
            Tolerance for acceptable deviation of each passed geometry
            coordinate  from that in the instance to still be considered
            matching. Default value is
            :data:`DEF.HESS_COORD_MATCH_TOL
            <opan.const.DEF.HESS_COORD_MATCH_TOL>`


        See :func:`opan.utils.check_geom <opan.utils.base.check_geom>` for
        details on return values and exceptions raised.

        """

        # Import(s)
        from .utils import check_geom as ucg

        # Wrapper call
        result = ucg(self.geom, self.atom_syms, coords, atoms, tol=tol)

        # Return result
        return result

    ## end def check_geom

## end class SuperOpanHess



class OrcaHess(SuperOpanHess):
    """ Container for HESS data generated by |orca|.

    Initialize by passing the path to the file to be loaded as the
    `path` keyword argument.

    Information contained includes the Hessian matrix, the number of atoms,
    the atomic symbols, the atomic weights, and the geometry, as reported in
    the .hess file.  See
    :ref:`Instance Variables <hess-orcahess-instancevars>`
    below for a full list.
    For variables marked *"required"* (those that should be found in all
    HESS files), a :exc:`~opan.error.HessError` is raised
    if the block is not found, whereas for all other variables
    a |None| value is stored.  For either type, if the data is malformed
    or invalid in some fashion, a :exc:`~opan.error.HessError` is raised with
    an appropriate typecode.

    Zero frequencies corresponding to translation/rotation are **not**
    excised from the frequencies list, normal modes, IR spectrum, Raman
    spectrum, etc.

    The precision of the geometry is less than that reported in an .xyz file,
    and thus should **not** be used for generation of subsequent computations.

    Output |units|:

         *  Hessian : Hartrees per Bohr-squared
            :math:`\\left(\\frac{\\mathrm{E_h}}{\\mathrm{B}^2}\\right)`

         *  Frequencies :  Wavenumbers
            :math:`\\left(\\frac{\\mathrm{cyc}}{\\mathrm{cm}}\\right)`

         *  IR intensities (:math:`T^2` values) :
            :math:`\\frac{\\mathrm{km}}{\\mathrm{mol}}`

         *  Raman activities :
            :math:`\\frac{\\mathring{\\mathrm{A}}^4}{\\mathrm u}`

         *  Dipole derivatives :
            :math:`\\mathrm{D_a\\over B} =
            \\mathrm{q\\, B\\over B} =
            \\mathrm{q}`

         *  Polarizability derivatives :
            :math:`\\mathrm{B^3\\over B} = \mathrm{B^2}`

         *  Eigenvalues of the mass-weighted Hessian :
            :math:`\\frac{\\mathrm{E_h}}{\\mathrm{u\\,B^2}}`

         *  Eigenvectors of the mass-weighted Hessian : (dimensionless)

    |

    **Methods**

    .. automethod:: _load

    |

    **Class Variables**

    .. class:: Pat

        |re.compile| patterns for data parsing.

        |

        .. attribute:: at_block

            Captures the entire block of atom IDs & weights, as
            well as the geometry data.

        .. attribute:: at_line

            Retrieves individual lines within the atom / geometry
            specification block.

        .. attribute:: dipder_block

            Captures the entire dipole derivatives block.

        .. attribute:: dipder_line

            Retrieves individual lines in the dipole derivatives block.

        .. attribute:: eigvals_block

            Captures the entire mass-weighted Hessian eigenvalues block.

        .. attribute:: eigvals_line

            Retrieves individual lines of the mass-weighted
            hessian eigenvalues block.

        .. attribute:: eigvecs_block

            Captures the entire set of mass-weighted Hessian
            eigenvectors blocks.

        .. attribute:: eigvecs_sec

            Extracts sections (individual blocks) of the mass-weighted
            Hessian eigenvectors blocks.

        .. attribute:: eigvecs_line

            Retrieves individual lines of a mass-weighted
            Hessian eigenvectors block.

        .. attribute:: energy

            Retrieves the energy reported in the .hess file.

        .. attribute:: freq_block

            Captures the entire vibrational frequencies block.

        .. attribute:: freq_line

            Retrieves individual lines of the frequencies block.

        .. attribute:: hess_block

            Captures the entire set of Hessian blocks.

        .. attribute:: hess_sec

            Extracts full-height, 3- or 6-column sections (individual
            blocks) of the set of Hessian blocks.

        .. attribute:: hess_line

            Retrieves individual lines within a Hessian block.

        .. attribute:: jobs_block

            Captures the entire job list block.

        .. attribute:: jobs_line

            Retrieves individual lines of the job list block.

        .. attribute:: ir_block

            Captures the entire IR spectrum block

        .. attribute:: ir_line

            Retrieves individual lines of the IR spectrum block

        .. attribute:: modes_block

            Captures the entire set of (weighted, column-normalized)
            normal modes blocks

        .. attribute:: modes_sec

            Extracts full-height, 3- or 6-column sections (individual
            blocks) of the set of normal modes blocks

        .. attribute:: modes_line

            Retrieves individual lines within a normal modes block

        .. attribute:: polder_block

            Captures the polarizability derivatives block

        .. attribute:: polder_line

            Retrieves individual lines in the polarizability
            derivatives block

        .. attribute:: raman_block

            Captures the Raman spectrum block

        .. attribute:: raman_line

            Retrieves individual lines in the Raman spectrum block

        .. attribute:: temp

            Retrieves the 'actual temperature' field (sometimes
            is spuriously zero)

    |

    .. _hess-orcahess-instancevars:

    **Instance Variables**

    .. attribute:: OrcaHess.atom_masses

        length-N |list| of |npfloat_| --
        List of atom masses as reported in the .hess file

    .. attribute:: OrcaHess.atom_syms

        length-N |list| of |str| *(required)* --
        List of uppercase atomic symbols

    .. attribute:: OrcaHess.dipders

        3N x 3 |npfloat_| --
        Matrix of dipole derivatives

    .. attribute:: OrcaHess.energy

        |float| *(required)* --
        Electronic energy reported in the Hessian file

    .. attribute:: OrcaHess.freqs

        length-3N |npfloat_| *(required)* --
        Vibrational frequencies in
        :math:`\\mathrm{cyc\\over cm}`,
        as reported in the Hessian file

    .. attribute:: OrcaHess.geom

        length-3N |npfloat_| *(required)* --
        Geometry vector

    .. attribute:: OrcaHess.hess

        3N x 3N |npfloat_| *(required)* --
        Cartesian Hessian matrix

    .. attribute:: OrcaHess.hess_path

        |str| --
        Complete path/filename from which the Hessian data was retrieved

    .. attribute:: OrcaHess.in_str

        |str| -- Complete contents of the imported .hess file

    .. attribute:: OrcaHess.ir_comps

        3N x 3 |npfloat_| --
        :math:`\\left(T_x, T_y, T_z\\right)`
        components of the transition dipole for each normal mode

    .. attribute:: OrcaHess.ir_mags

        length-3N |npfloat_| --
        :math:`T^2` values (squared-magnitudes) of the
        transition dipole for each mode

    .. attribute:: OrcaHess.joblist

        N x 3 |bool| --
        Completion status for each displacement in calculation of the Hessian

    .. attribute:: OrcaHess.modes

        3N x 3N |npfloat_| *(required)* --
        Rotation- and translation-purified, mass-weighted
        vibrational normal modes,
        with each mode (column vector) separately normalized by |orca|.

    .. attribute:: OrcaHess.mwh_eigvals

        length-3N |npfloat_| --
        Eigenvalues of the mass-weighted Hessian, in
        :math:`\\mathrm{E_h \\over u\\,B^2}`

    .. attribute:: OrcaHess.mwh_eigvecs

        3N x 3N |npfloat_| --
        Eigenvectors of the mass-weighted Hessian, as column vectors: the
        eigenvector of eigenvalue :math:`i` would be retrieved with
        :samp:`mwh_eigvecs[:,{i}]`

    .. attribute:: OrcaHess.num_ats

        |int| -- Number of atoms in the system

    .. attribute:: OrcaHess.polders

        3N x 6 |npfloat_| --
        Matrix of Cartesian polarizability derivatives

    .. attribute:: OrcaHess.raman_acts

        length-3N |npfloat_| --
        Vector of Raman activities

    .. attribute:: OrcaHess.raman_depols

        length-3N |npfloat_| --
        Vector of Raman depolarization factors

    .. attribute:: OrcaHess.temp

        |float| *(required)* --
        "Actual temperature" reported in the .hess file. Occasionally
        stored by |orca| as a meaningless zero value instead
        of the temperature used.

    """

    # Imports

    # Various class-level RegEx patterns, wrapped in a class
    class Pat(object):
        """ Various class-level RegEx patterns, wrapped in a class. """

        # Imports
        import re as _re

        # Various class-level Regex patterns.

        # Atoms list, including atomic weights. Assumes no scientific notation
        #  will be used in the coordinates.
        at_block = _re.compile("""
        \\#.*\\n                    # Line prior to $atoms is a blank comment
        \\$atoms.*\\n                   # $atoms indicator for start of block
        (?P<num>[0-9]+).*               # Storing the number of atoms
        (?P<block>\\n                   # Store the whole chunk of coordinates
            (                           # Open group for atom lines def
                [ \\t]+([a-z]+|[0-9]+)  # Whitespace and atomic number/element
                [ \\t]+[0-9.]+          # Whitespace and atomic mass
                (                       # Open group for coordinates
                    [ \\t]+[0-9.-]+     # Whitespace and one coordinate
                ){3}                    # Three coordinates for each line
                .*\\n                   # Whatever to end of line
            )+                          # Some number of coordinate lines
        )                               # Close the "block" group
        """, _re.I | _re.M | _re.X)

        # Pulling the individual atoms, weights, and coordinates from each atom line
        at_line = _re.compile("""
        ^[ \\t]+                        # Whitespace to start the line
        (?P<el>([a-z]+|[0-9]+))         # Atomic number or element symbol
        [ \\t]+(?P<mass>[0-9.]+)        # Whitespace and atomic mass
        [ \\t]+(?P<c1>[0-9.-]+)         # Whitespace and first coordinate
        [ \\t]+(?P<c2>[0-9.-]+)         # Whitespace and second coordinate
        [ \\t]+(?P<c3>[0-9.-]+)         # Whitespace and third coordinate
        """, _re.I | _re.M | _re.X)

        # Entire Hessian data block
        hess_block = _re.compile("""
        \\$hessian.*\\n             # Marker for Hessian block
        (?P<dim>[0-9]+).*\\n        # Dimensionality of Hessian (3N x 3N)
        (?P<block>                  # Group for the subsequent block of lines
            (                       # Group for single line definition
                ([ \\t]+[0-9.-]+)+  # Some number of whitespace-separated nums
                .*\\n               # Plus whatever to end of line
            )+                      # Whatever number of single lines
        )                           # Enclose the whole batch of lines
        """, _re.I | _re.X)

        # Sections of the Hessian data block
        hess_sec = _re.compile("""
        ([ \\t]+[0-9]+)+[ \\t]*\\n  # Column header line
        (                           # Open the group for the sub-block lines
            [ \\t]+[0-9]+           # Row header
            (                       # Open the group defining a single element
                [ \\t]+[-]?         # Whitespace and optional hyphen
                [0-9]+\\.[0-9]+     # One or more digits, decimal, more digits
            )+                      # Some number of sub-columns
            [ \\t]*\\n              # Whitespace to EOL
        )+                          # Some number of suitable lines
        """, _re.I | _re.X)

        # Pulling Hessian lines from the sections, with elements in groups
        hess_line = _re.compile("""
        ^[ \\t]*                        # Optional whitespace to start each line
        (?P<row>[0-9]+)                         # Row header
        [ \\t]+(?P<e0>[0-9-]+\\.[0-9]+)         # 1st element
        [ \\t]+(?P<e1>[0-9-]+\\.[0-9]+)         # 2nd element
        [ \\t]+(?P<e2>[0-9-]+\\.[0-9]+)         # 3rd element
        ([ \\t]+(?P<e3>[0-9-]+\\.[0-9]+))?      # 4th element (possibly absent)
        ([ \\t]+(?P<e4>[0-9-]+\\.[0-9]+))?      # 5th element (possibly absent)
        ([ \\t]+(?P<e5>[0-9-]+\\.[0-9]+))?      # 6th element (possibly absent)
        .*$                             # Whatever to end of line
        """, _re.I | _re.M | _re.X)

        # Reported energy
        energy = _re.compile("""
        \\$act_energy[ ]*\\n                    # Label for the block
        [ ]+(?P<en>[0-9.-]+)[ ]*\\n             # Energy value
        """, _re.I | _re.X)

        # Frequencies block
        freq_block = _re.compile("""
        \\$vibrational_frequencies[ ]*\\n       # Block label
        [ ]*(?P<num>[0-9]+)[ ]*\\n              #--> Number of frequencies
        (?P<block>                              # Open capture group for block
            ([ ]+[0-9]+[ ]+[0-9.-]+[ ]*\\n)+    #--> Block contents
        )                                       # Close capture group
        """, _re.I | _re.X)

        freq_line = _re.compile("""
        ^[ ]+(?P<id>[0-9]+)                     #--> ID for the frequency
        [ ]+(?P<freq>[0-9.-]+)                  #--> Freq value in cm**-1
        [ ]*$                                   # To EOL
        """, _re.I | _re.M | _re.X)


        # Entire modes data block
        modes_block = _re.compile("""
        \\$normal_modes.*\\n        # Marker for modes block
        (?P<dim>[0-9]+)[ ]+         # Dimensionality of modes block (3N x 3N)
        (?P<dim2>[0-9]+).*\\n       #  (Second dimension value)
        (?P<block>                  # Group for the subsequent block of lines
            (                       # Group for single line definition
                ([ ]+[0-9.-]+)+     # Some number of whitespace-separated nums
                .*\\n               # Plus whatever to end of line
            )+                      # Whatever number of single lines
        )                           # Enclose the whole batch of lines
        """, _re.I | _re.X)

        # Sections of the modes data block
        modes_sec = _re.compile("""
        ([ ]+[0-9]+)+[ ]*\\n        # Column header line
        (                           # Open the group for the sub-block lines
            [ ]+[0-9]+              # Row header
            (                       # Open the group defining a single element
                [ ]+[-]?            # Whitespace and optional hyphen
                [0-9]+\\.[0-9]+     # One or more digits, decimal, more digits
            )+                      # Some number of sub-columns
            [ ]*\\n                 # Whitespace to EOL
        )+                          # Some number of suitable lines
        """, _re.I | _re.X)

        # Pulling modes lines from the sections, with elements in groups
        #  THE USE of the '[0-9-]+\\.[0-9]+' construction here, and in its variants
        #  above/below, guarantees a floating-point value is found. Otherwise, the
        #  Regex retrieves on into subsequent sections because the header rows
        #  parse just fine for a '[ ]+[0-9.-]' pattern.
        modes_line = _re.compile("""
        ^[ ]*                           # Optional whitespace to start each line
        (?P<row>[0-9]+)                         # Row header
        [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            # 1st element
        [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            # 2nd element
        [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            # 3rd element
        ([ ]+(?P<e3>[0-9-]+\\.[0-9]+))?         # 4th element (possibly absent)
        ([ ]+(?P<e4>[0-9-]+\\.[0-9]+))?         # 5th element (possibly absent)
        ([ ]+(?P<e5>[0-9-]+\\.[0-9]+))?         # 6th element (possibly absent)
        .*$                             # Whatever to end of line
        """, _re.I | _re.M | _re.X)

        # "Actual temperature"
        temp = _re.compile("""
        \\$actual_temperature[ ]*\\n            # Marker for value
        [ ]*(?P<temp>[0-9.]+)[ ]*\\n            # Value
        """, _re.I | _re.X)

        # Dipole derivatives block
        dipder_block = _re.compile("""
        \\$dipole_derivatives[ ]*\\n            # Marker for block
        (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
        (?P<block>                              #--> Catches entire block
            (([ ]+[0-9.-]+)+[ ]*\\n)+           # Rows of data
        )                                       # End block capture
        """, _re.I | _re.X)

        # Dipole derivatives individual line
        dipder_line = _re.compile("""
        ^                                       # Line start
        [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            # 1st element
        [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            # 2nd element
        [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            # 3rd element
        [ ]*$
        """, _re.I | _re.M | _re.X)

        # IR spectrum block
        ir_block = _re.compile("""
        \\$ir_spectrum[ ]*\\n                   # Marker for block
        (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
        (?P<block>                              #--> Catch entire block
            (([ ]+[0-9.-]+)+[ ]*\\n)+           # Rows of data
        )                                       # End block catch
        """, _re.I | _re.X)

        ir_line = _re.compile("""
        ^                                       # Line start
        [ ]+(?P<freq>[0-9-]+\\.[0-9]+)          #--> Frequency
        [ ]+(?P<mag>[0-9]+\\.[0-9]+)            #--> Transition dipole sq. mag
        [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            #--> 1st element (TX)
        [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            #--> 2nd element (TY)
        [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            #--> 3rd element (TZ)
        """, _re.I | _re.M | _re.X)

        # Polarizability derivatives block
        polder_block = _re.compile("""
        \\$polarizability_derivatives[ ]*\\n    # Block marker
        (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
        (?P<block>                              #--> Catch entire block
            (([ ]+[0-9.-]+)+[ ]*\\n)+           # Rows of data
        )                                       # End block catch
        """, _re.I | _re.X)

        polder_line = _re.compile("""
        ^                                       # Start of line
        [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            # 1st element
        [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            # 2nd element
        [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            # 3rd element
        [ ]+(?P<e3>[0-9-]+\\.[0-9]+)            # 4th element
        [ ]+(?P<e4>[0-9-]+\\.[0-9]+)            # 5th element
        [ ]+(?P<e5>[0-9-]+\\.[0-9]+)            # 6th element
        """, _re.I | _re.M | _re.X)


        # Raman spectrum block
        raman_block = _re.compile("""
        \\$raman_spectrum[ ]*\\n                # Block marker
        (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
        (?P<block>                              #--> Catch entire block
            (([ ]+[0-9.-]+)+[ ]*\\n)+           # Rows of data
        )                                       # End block catch
        """, _re.I | _re.X)

        raman_line = _re.compile("""
        ^                                       # Start of line
        [ ]+(?P<freq>[0-9-]+\\.[0-9]+)          #--> Frequency of mode
        [ ]+(?P<act>[0-9]+\\.[0-9]+)            #--> Raman activity
        [ ]+(?P<depol>[0-9]+\\.[0-9]+)          #--> Depolarization factor
        """, _re.I | _re.M | _re.X)


        # Job list block
        jobs_block = _re.compile("""
        \\$job_list[ ]*\\n                      # Block marker
        (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
        (?P<block>                              #--> Catch entire block
            (([ ]+[0-9]+)+[ ]*\\n)+             # Rows of data
        )                                       # End block catch
        """, _re.I | _re.X)

        jobs_line = _re.compile("""
        ^                                       # Start of line
        [ ]+(?P<at>[0-9]+)                      #--> Atom index
        [ ]+(?P<e0>[01])                        # 1st element
        [ ]+(?P<e1>[01])                        # 2nd element
        [ ]+(?P<e2>[01])                        # 3rd element
        """, _re.I | _re.M | _re.X)

        # Eigenvalues of mass-weighted Hessian
        eigvals_block = _re.compile("""
        \\$eigenvalues_mass_weighted_hessian[ ]*\\n         # Block marker
        (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
        (?P<block>                              #--> Catch entire block
            (([ ]+[0-9.-]+)+[ ]*\\n)+           # Rows of data
        )                                       # End block catch
        """, _re.I | _re.X)

        eigvals_line = _re.compile("""
        ^                                       # Start of line
        [ ]+(?P<mode>[0-9]+)                    #--> Mode index
        [ ]+(?P<eig>[0-9.-]+)                   #--> Eigenvalue
        """, _re.I | _re.M | _re.X)


        #=== Mass-weighted Hessian eigenvectors ===#
        # Entire eigenvectors data block
        eigvecs_block = _re.compile("""
        \\$eigenvectors_mass_weighted_hessian.*\\n  # Marker for modes block
        (?P<dim>[0-9]+)[ ]+         # Dimensionality of modes block (3N x 3N)
        (?P<dim2>[0-9]+)[ ]*\\n     #  (Second dimension value)
        (?P<block>                  # Group for the subsequent block of lines
            (                       # Group for single line definition
                ([ ]+[0-9.-]+)+     # Some number of whitespace-separated nums
                .*\\n               # Plus whatever to end of line
            )+                      # Whatever number of single lines
        )                           # Enclose the whole batch of lines
        """, _re.I | _re.X)

        # Sections of the eigenvectors data block
        eigvecs_sec = _re.compile("""
        ([ ]+[0-9]+)+[ ]*\\n        # Column header line
        (                           # Open the group for the sub-block lines
            [ ]+[0-9]+              # Row header
            (                           # Open the group defining a single element
                [ ]+[-]?                # Whitespace and optional hyphen
                [0-9]+\\.[0-9]+         # One or more digits, decimal, more digits
            )+                          # Some number of sub-columns
            [ ]*\\n                     # Whitespace to EOL
        )+                              # Some number of suitable lines
        """, _re.I | _re.X)

        # Pulling modes lines from the sections, with elements in groups
        #  THE USE of the '[0-9-]+\\.[0-9]+' construction here, and in its variants
        #  above/below, guarantees a floating-point value is found. Otherwise, the
        #  Regex retrieves on into subsequent sections because the header rows
        #  parse just fine for a '[ ]+[0-9.-]' pattern.
        eigvecs_line = _re.compile("""
        ^[ ]*                           # Optional whitespace to start each line
        (?P<row>[0-9]+)                         # Row header
        [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            # 1st element
        [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            # 2nd element
        [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            # 3rd element
        ([ ]+(?P<e3>[0-9-]+\\.[0-9]+))?         # 4th element (possibly absent)
        ([ ]+(?P<e4>[0-9-]+\\.[0-9]+))?         # 5th element (possibly absent)
        ([ ]+(?P<e5>[0-9-]+\\.[0-9]+))?         # 6th element (possibly absent)
        .*$                             # Whatever to end of line
        """, _re.I | _re.M | _re.X)

    ## end class Pat


    def _load(self, **kwargs):
        """ Initialize OrcaHess Hessian object from .hess file

        Searches indicated file for data blocks within the .hess file.  The
        geometry, Hessian block, frequencies, and normal modes must be present
        and will be retrieved; other blocks will be imported if present, or
        ignored if absent (stored as |None| in the resulting object).
        If malformed/inaccurate data is found in any
        block that is present, some flavor of :exc:`~opan.error.HessError`
        will be raised.

        Parameters
        ----------
        path
            |str| --
            `kwargs` parameter specifying complete path to the
            .hess file to be read.

        Raises
        ------
        ~opan.error.HessError
            (various typecodes) If indicated Hessian file
            is malformed in some fashion.

        ~excoptions.KeyError
            If invalid atomic symbol appears in .hess file.

        ~exceptions.IOError
            If the indicated file does not exist or cannot be read.


        """

        # Local method(s)
        def parse_multiblock(hesstext, p_block, p_sec, p_line, num_ats,
                                                            blockname, tc):
            """ Helper function for importing blocks with multiple sections.

            Parsing of data spanning multiple sections of columns is somewhat
            involved. This function encapsulates the process for cleaner
            code.  The structure depends critically on several formatting
            features of |orca| .hess files.

            Note the search groups that must be present in the `p_block` and
            `p_line` Regex patterns.

            Parameters
            ----------
            hesstest
                |str| --
                Complete text of the .hess file

            p_block
                |re.compile| pattern --
                Retrieves the **entirety** of the relevant block
                Required groups:
                    dim     : overall dimension of the data block
                    block   : contents of the block, including column headers

            p_sec
                |re.compile| pattern --
                Retrieves each section of data columns

            p_line
                |re.compile| pattern --
                Retrieves individual lines of a section
                Required groups:
                    row     : Row index into the final data block
                    e#      : Column index into the data section (# in range(6))

            num_ats
                |int| --
                Number of atoms in the geometry

            blockname
                |str| --
                Brief text description of the block being imported, if needed
                for error reporting purposes

            tc
                :class:`~opan.error.HessError` typecode --
                Type of error to be thrown, if required

            Returns
            -------
            workmtx
                3N x 3N |npfloat_| --
                Assembled array for storage

            """

            # Pull the block
            m_block = p_block.search(hesstext)

            # Confirm the anticipated matrix size matches 3*num_ats
            if not int(m_block.group("dim")) == 3*num_ats:
                raise HessError(tc, blockname + " dimension " +
                        "specification mismatched with geometry",
                        "HESS File: " + hess_path)
            ## end if

            # Initialize the working matrix
            workmtx = np.zeros((3*num_ats, 3*num_ats), dtype=np.float_)

            # Initialize the column offset for populating the matrix
            col_offset = 0

            # Initialize the counter for the number of row sections imported
            rows_counter = 0

            # Loop through the subsections of the matrix
            for m_sec in p_sec.finditer(m_block.group("block")):
                # Loop through each entry
                for m_line in p_line.finditer(m_sec.group(0)):
                    # Store the row; expect base zero so no adjustment needed.
                    rowval = scast(m_line.group("row"), np.int_)

                    # Loop to fill the row
                    for i in range(6):
                        # Store the element
                        val = scast(m_line.group('e{0}'.format(i)), np.float_)

                        # Only store to matrix if a value actually retrieved.
                        #  This protects against the final three-column section
                        #  in blocks for systems with an odd number of atoms.
                        if not np.isnan(val):
                            workmtx[rowval, col_offset + i] = val
                        ## end if
                    ## next i

                    # Increment the row-read counter
                    rows_counter += 1

                ## next m_hess_line

                # Last thing is to increment the offset by six. Don't have to
                #  worry about offsetting for a section only three columns wide
                #  because that SHOULD only ever occur at the last section
                col_offset += 6

            ## next m_hess_sec

            # Check to ensure that the column offset is high enough to have
            #  fully populated the matrix based on the number of atoms. This is
            #  a check against malformation of the HESS file, resulting in a
            #  reduced number of sections being retrieved.
            if not col_offset >= (3 * num_ats):
                raise HessError(tc, "Insufficient number of " + blockname +
                        " sections found",
                        "HESS File: " + hess_path)
            ## end if

            # Additional cross-check on number of rows imported
            if (rows_counter % (3*num_ats)) != 0:
                # Not the right number of rows; complain
                raise HessError(tc,
                            blockname + " row count mismatch",
                            "HESS File: " + hess_path)
            ## end if

            # Return the working matrix
            return workmtx

        ## end def parse_multiblock

        # Imports
        import re, numpy as np
        from .error import HessError
        from .utils import safe_cast as scast
        from .const import atom_num, atom_sym, PRM, DEF

        # Check if instantiated; complain if so
        if 'hess_path' in dir(self):
            raise HessError(HessError.OVERWRITE,
                    "Cannot overwrite contents of existing OrcaHess", "")
        ## end if

        # Retrieve the file target and store
        hess_path = kwargs['path']
        self.hess_path = hess_path

        # Open file, read contents, close stream
        with open(hess_path,'rU') as in_fl:
            self.in_str = in_fl.read()
        ## end with

        # Store the source string
        srcstr = "HESS File: {0}".format(hess_path)

        # Check to ensure all required data blocks are found
        if not OrcaHess.Pat.at_block.search(self.in_str):
            raise HessError(HessError.AT_BLOCK,
                    "Atom specification block not found", srcstr)
        ## end if
        if not OrcaHess.Pat.hess_block.search(self.in_str):
            raise HessError(HessError.HESS_BLOCK,
                    "Hessian block not found", srcstr)
        ## end if
        if not OrcaHess.Pat.freq_block.search(self.in_str):
            raise HessError(HessError.FREQ_BLOCK,
                    "Frequencies block (cm**-1 units) not found", srcstr)
        ## end if
        if not OrcaHess.Pat.modes_block.search(self.in_str):
            raise HessError(HessError.MODES_BLOCK,
                    "Normal modes block not found", srcstr)
        ## end if
        if not OrcaHess.Pat.energy.search(self.in_str):
            raise HessError(HessError.ENERGY,
                    "Energy value not found", srcstr)
        ## end if
        if not OrcaHess.Pat.temp.search(self.in_str):
            raise HessError(HessError.TEMP,
                    "'Actual temperature' value not found", srcstr)
        ## end if


        #=== Geometry spec block ===#
        # Store the block
        m_work = self.Pat.at_block.search(self.in_str)

        # Bring in the number of atoms
        self.num_ats = np.int_(m_work.group("num"))

        # Initialize the vectors of atomic symbols, masses, and the
        #  geometry
        self.atom_syms = []
        self.atom_masses = []
        self.geom = []

        # Parse the list of atoms
        for m in OrcaHess.Pat.at_line.finditer(m_work.group("block")):
            # Parse the element symbol or atomic number
            try:
                # See if it casts as an int
                num = scast(m.group('el'), np.int_)
            except (ValueError, TypeError):
                # Nope, must be letters. Check if valid by poking it into the
                #  dict. If not valid, should raise another error.
                num = atom_num[m.group('el').upper()]
            ## end try

            # If this point is reached, either it's cast ok to an int,
            #  or the string has been passed through the atom_num dict and
            #  is a valid element symbol. SO, convert back to atom_sym and
            #  store - this will also check for an invalid int entry
            self.atom_syms.append(atom_sym[num])

            # Parse the atomic weight; should be a simple append?
            self.atom_masses.append(scast(m.group("mass"), np.float_))

            # Build geometry as a 1-D list
            self.geom.extend([scast(m.group(gp), np.float_) for
                                                gp in ['c1', 'c2', 'c3']])
        ## next m

        # Double-check that the number of atoms retrieved matches the
        #  number indicated in the HESS file. Checks of atomic masses
        #  and geometry length are presumably redundant.
        if not len(self.atom_syms) == self.num_ats:
            raise HessError(HessError.AT_BLOCK, "Atomic symbol dimension " +
                    "mismatch with HESS atom specification", srcstr)
        ## end if

        # Convert geometry to numpy storage form
        self.geom = np.array(self.geom, dtype=np.float_)

        #=== Hessian and modes ===#
        # Pull the Hessian
        self.hess = parse_multiblock(self.in_str, self.Pat.hess_block,
                self.Pat.hess_sec, self.Pat.hess_line, self.num_ats,
                "Hessian", HessError.HESS_BLOCK)

        # Pull the modes
        self.modes = parse_multiblock(self.in_str, self.Pat.modes_block,
                self.Pat.modes_sec, self.Pat.modes_line, self.num_ats,
                "modes", HessError.MODES_BLOCK)

        # Extra check of 'dim' vs 'dim2' on modes
        m_work = self.Pat.modes_block.search(self.in_str)
        if scast(m_work.group("dim"), np.int_) != \
                            scast(m_work.group("dim2"), np.int_):
            raise HessError(HessError.MODES_BLOCK,
                    "Normal mode block dimension specification mismatch",
                    srcstr)
        ## end if

        #=== Frequencies ===#
        # Pull the block
        m_work = self.Pat.freq_block.search(self.in_str)

        # Check that number of frequencies indicated in the block matches
        #  that expected from the number of atoms
        if 3*self.num_ats != scast(m_work.group("num"), np.int_):
            raise HessError(HessError.FREQ_BLOCK,
                    "Count in frequencies block != 3 * number of atoms",
                    srcstr)
        ## end if

        # Retrieve the frequencies
        self.freqs = np.array(
                [scast(m.group("freq"), np.float_) for m in
                self.Pat.freq_line.finditer(m_work.group("block"))])

        # Proofread for proper size
        if not self.freqs.shape[0] == 3*self.num_ats:
            raise HessError(HessError.FREQ_BLOCK,
                    "Number of frequencies != 3 * number of atoms", srcstr)
        ## end if


        #=== Pull the single values that should always be present ===#
        # Store the reported energy
        self.energy = scast(self.Pat.energy.search(self.in_str).group("en"),
                                                                    np.float_)

        # Store the reported 'actual temperature'
        self.temp = scast(self.Pat.temp.search(self.in_str).group("temp"),
                                                                    np.float_)

        #=== Dipole derivatives ===#
        # Check if block found. Store None if not; otherwise import
        m_work = self.Pat.dipder_block.search(self.in_str)
        if m_work is None:
            self.dipders = None
        else:
            # Check that number of derivatives rows indicated in the block
            #  matches that expected from the number of atoms
            if 3*self.num_ats != scast(m_work.group("dim"), np.int_):
                raise HessError(HessError.DIPDER_BLOCK,
                        "Count in dipole derivatives block != 3 * # of atoms",
                        srcstr)
            ## end if

            # Retrieve the derivatives
            self.dipders = np.array(
                    [[scast(m.group("e{0}".format(i)), np.float_)
                        for i in range(3)]
                        for m in self.Pat.dipder_line.finditer(
                                                    m_work.group("block")) ])

            # Proofread for proper size. Don't have to proofread the width of
            #  three, since any row not containing three numerical values will
            #  result in the block getting truncated, per the Regex constrution.
            if not self.dipders.shape[0] == 3*self.num_ats:
                raise HessError(HessError.DIPDER_BLOCK,
                        "Number of dipole derivative rows != \
                        3 * number of atoms", srcstr)
            ## end if

            # If max-absolute element is too big, overwrite with None
            if np.max(np.abs(self.dipders)) > PRM.MAX_SANE_DIPDER:
                self.dipders = None
            ## end if
        ## end if


        #=== IR Spectrum ===#
        # If dipole derivs absent or munged, or if block missing, then skip
        m_work = self.Pat.ir_block.search(self.in_str)
        if self.dipders is None or m_work is None:
            self.ir_comps = None
            self.ir_mags = None
        else:
            # Complain if number of stated modes mismatches expectation
            if 3*self.num_ats != np.int_(m_work.group("dim")):
                raise HessError(HessError.IR_BLOCK,
                        "Count in IR spectrum block != 3 * # of atoms",
                        srcstr)
            ## end if

            # Pull the blocks
            self.ir_comps = np.array(
                    [[scast(m.group("e{0}".format(i)), np.float_)
                        for i in range(3)]
                        for m in self.Pat.ir_line.finditer(
                                                    m_work.group("block")) ])
            self.ir_mags = np.array(
                    [scast(m.group("mag"), np.float_)
                        for m in self.Pat.ir_line.finditer(
                                                    m_work.group("block")) ])

            # Confirm length of ir_mags conforms. Shouldn't need to check both,
            #  since they both rely equally on Pat.ir_line.finditer.
            if 3*self.num_ats != self.ir_mags.shape[0]:
                raise HessError(HessError.IR_BLOCK,
                        "Number of IR spectrum rows != \
                        3 * number of atoms", srcstr)
            ## end if

            # Confirm match of all frequencies with those reported separately
            if not np.allclose(
                    self.freqs,
                    np.array([np.float_(m.group('freq'))
                                for m in self.Pat.ir_line.finditer(
                                m_work.group("block") )]),
                    rtol=0,
                    atol=DEF.HESS_IR_MATCH_TOL):
                raise HessError(HessError.IR_BLOCK,
                        "Frequency mismatch between freq and IR blocks",
                        srcstr)
            ## end if
        ## end if


        #=== Polarizability Derivatives ===#
        # If block is missing, skip it
        m_work = self.Pat.polder_block.search(self.in_str)
        if m_work is None:
            self.polders = None
        else:
            # Check that number of derivatives rows indicated in the block
            #  matches that expected from the number of atoms
            if 3*self.num_ats != np.int_(m_work.group("dim")):
                raise HessError(HessError.POLDER_BLOCK,
                        "Count in polarizability derivatives block \
                        != 3 * # of atoms", srcstr)
            ## end if

            # Retrieve the derivatives
            self.polders = np.array(
                    [[scast(m.group("e{0}".format(i)), np.float_)
                        for i in range(6)]
                        for m in self.Pat.polder_line.finditer(
                                                    m_work.group("block")) ])

            # Proofread for proper size. Don't have to proofread the width of
            #  six, since any row not containing six numerical values will
            #  result in the block getting truncated.
            if not self.polders.shape[0] == 3*self.num_ats:
                raise HessError(HessError.POLDER_BLOCK,
                        "Number of polarizability derivative rows != \
                        3 * number of atoms", srcstr)
            ## end if
        ## end if


        #=== Raman Spectrum ===#
        # If polarizability derivs absent or munged, or if block missing,
        #  then skip
        m_work = self.Pat.raman_block.search(self.in_str)
        if self.polders is None or m_work is None:
            self.raman_acts = None
            self.raman_depols = None
        else:
            # Complain if number of stated modes mismatches expectation
            if 3*self.num_ats != np.int_(m_work.group("dim")):
                raise HessError(HessError.RAMAN_BLOCK,
                        "Count in Raman spectrum block != 3 * # of atoms",
                        srcstr)
            ## end if

            # Pull the blocks
            self.raman_acts = np.array(
                    [np.float_(m.group("act"))
                        for m in self.Pat.raman_line.finditer(
                                                    m_work.group("block")) ])
            self.raman_depols = np.array(
                    [np.float_(m.group("depol"))
                        for m in self.Pat.raman_line.finditer(
                                                    m_work.group("block")) ])

            # Confirm length of raman_acts conforms. Shouldn't need to check
            #  both, since they both rely equally on Pat.raman_line.finditer.
            if 3*self.num_ats != self.raman_acts.shape[0]:
                raise HessError(HessError.RAMAN_BLOCK,
                        "Number of Raman spectrum rows != \
                        3 * number of atoms", srcstr)
            ## end if

            # Confirm match of all frequencies with those reported separately
            if not np.allclose(
                    self.freqs,
                    np.array([np.float_(m.group('freq'))
                                for m in self.Pat.raman_line.finditer(
                                m_work.group("block")) ]),
                    rtol=0,
                    atol=DEF.HESS_IR_MATCH_TOL):
                raise HessError(HessError.RAMAN_BLOCK,
                        "Frequency mismatch between freq and Raman blocks",
                        srcstr)
            ## end if
        ## end if


        #=== Job list ===#
        # Check if block found. Store None if not; otherwise import
        m_work = self.Pat.jobs_block.search(self.in_str)
        if m_work is None:
            self.joblist = None
        else:
            # Check that number of joblist rows indicated in the block
            #  matches that expected from the number of atoms
            if 3*self.num_ats != np.int_(m_work.group("dim")):
                raise HessError(HessError.JOB_BLOCK,
                        "Count in job list block != 3 * # of atoms", srcstr)
            ## end if

            # Retrieve the job list
            self.joblist = np.array(
                    [[scast(m.group("e{0}".format(i)), np.float_)
                        for i in range(3)]
                        for m in self.Pat.jobs_line.finditer(
                                                    m_work.group("block")) ])

            # Proofread for proper size. Don't have to proofread the width of
            #  three, since any row not containing three numerical values will
            #  result in the block getting truncated.
            if not self.joblist.shape[0] == self.num_ats:
                raise HessError(HessError.JOB_BLOCK,
                        "Number of job list rows != number of atoms", srcstr)
            ## end if

            # Convert to boolean
            self.joblist = np.equal(self.joblist, np.ones_like(self.joblist))

        ## end if


        #=== Mass-weighted Hessian -- Eigenvalues ===#
        # Pull the block; continue only if block is present
        m_work = self.Pat.eigvals_block.search(self.in_str)
        if m_work is None:
            self.mwh_eigvals = None
        else:
            # Check that number of eigenvalues indicated in the block matches
            #  that expected from the number of atoms
            if 3*self.num_ats != np.int_(m_work.group("dim")):
                raise HessError(HessError.EIGVAL_BLOCK,
                        "Count in MWH eigenvalues block != 3 * number of atoms",
                        srcstr)
            ## end if

            # Retrieve the eigenvalues
            self.mwh_eigvals = np.array(
                    [np.float_(m.group("eig")) for m in
                    self.Pat.eigvals_line.finditer(m_work.group("block"))])

            # Proofread for proper size
            if not self.mwh_eigvals.shape[0] == 3*self.num_ats:
                raise HessError(HessError.EIGVAL_BLOCK,
                        "Number of MWH eigenvalues != 3 * number of atoms",
                        srcstr)
            ## end if

        ## end if


        #=== Mass-weighted Hessian -- Eigenvalues ===#
        # See if the block is there; import if so
        m_work = self.Pat.eigvecs_block.search(self.in_str)
        if m_work is None:
            self.mwh_eigvecs = None
        else:
            # Pull the eigenvectors
            self.mwh_eigvecs = \
                    parse_multiblock(self.in_str, self.Pat.eigvecs_block,
                    self.Pat.eigvecs_sec, self.Pat.eigvecs_line, self.num_ats,
                    "MWH eigenvectors", HessError.EIGVEC_BLOCK)

            # Extra check of 'dim' vs 'dim2' on modes
            if scast(m_work.group("dim"), np.int_) != \
                                scast(m_work.group("dim2"), np.int_):
                raise HessError(HessError.EIGVEC_BLOCK,
                        "MWH eigenvectors dimension specification mismatch",
                        srcstr)
            ## end if
        ## end if

    ## end def __init__

## end class OrcaHess


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")


