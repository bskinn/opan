#-------------------------------------------------------------------------------
# Name:        grad
# Purpose:     Encapsulation of parsing/handling of data from gradient files.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     28 Oct 2014
# Copyright:   (c) Brian Skinn 2015
# License:     The MIT License; see "license.txt" for full license terms
#                   and contributor agreement.
#
#       This file is part of opan (OpenAnharmonic), a system for automated
#       computation of anharmonic properties of molecular systems via wrapper
#       calls to computational/quantum chemical software packages.
#
#       http://www.github.com/bskinn/opan
#
#-------------------------------------------------------------------------------

""" Module implementing imports of gradient data from external computations.

The abstract superclass :class:`SuperOpanGrad` defines a common initializer
and common method(s) that are to be used by subclasses designed to import
gradient data from external computational packages.  The import for each
external software should have its own subclass.

Units of the gradient as stored should always be Hartrees per Bohr
:math:`\\left(\\frac{\\mathrm{E_h}}{\\mathrm B}\\right)`.

Implemented subclasses are:

:class:`OrcaEngrad` -- Imports '.engrad' files from |orca|

|

**Superclass**

.. autoclass:: SuperOpanGrad

|

**Subclasses**

.. autoclass:: OrcaEngrad(path='...')
    :members:

"""


# Imports


# Debug constant
_DEBUG = False


#TODO: Make abstract ENGRAD superclass and put check_geom into it; then make
#  all specific grad classes inherit from it.

class SuperOpanGrad(object):
    """ Abstract superclass of gradient import classes.

    Performs the following actions:

    1.  Ensures that the abstract superclass is not being instantiated,
        but instead a subclass.
    2.  Calls the ``_load()`` method on the subclass, passing all `kwargs`
        through unmodified.
    3.  Typechecks the ``self.gradient``, ``self.geom``, and
        ``self.atom_syms`` required members for existence, proper data type,
        and properly matched lengths.

    |

    **Methods**

    .. automethod:: check_geom(coords, atoms[, tol])

    """

    # Imports
    from .const import DEF as _DEF

    def __init__(self, **kwargs):
        """ Conformity-checking initalizer for gradient data loader classes.
        """

        # Imports
        from .error import GradError as GErr
        from .utils import assert_npfloatarray as a_npfa
        from .const import atomNum
        import numpy as np

        # Check for abstract base class
        if type(self) == SuperOpanGrad:
            raise(NotImplementedError("SuperOpanGrad base class is abstract"))
        ## end if

        # Call the subclass _load method, passing in the keyword arguments
        #  wholesale
        self._load(**kwargs)

        # Define common error source string
        srcstr = "{0} with args {1}".format(self.__class__, str(kwargs))

        # All proofreading here is excluded from coverage because
        #  all of these should be handled at the subclass level for ORCA.

        # Proofread types for gradient and geometry
        a_npfa(self, 'gradient', 'gradient', GErr, GErr.badgrad, srcstr)
        a_npfa(self, 'geom', 'geometry', GErr, GErr.badgeom, srcstr)

        # Gradient and geometry shape
        if (len(self.gradient.shape) != 1 or
                                    self.gradient.shape[0] % 3 != 0):
            raise(GErr(GErr.badgrad, # pragma: no cover
                        "Gradient is not a length-3N vector", srcstr))
        ## end if
        if self.gradient.shape != self.geom.shape:
            raise(GErr(GErr.badgeom, # pragma: no cover
                        "Geometry shape does not match gradient shape", srcstr))

        # Ensure atomic symbols length matches & they're all valid
        if not hasattr(self, 'atom_syms'): # pragma: no cover
            raise(GErr(GErr.badatom, "Atoms list not found", srcstr))
        ## end if
        if 3*len(self.atom_syms) != self.gradient.shape[0]: # pragma: no cover
            raise(GErr(GErr.badatom, "Atoms list is not length-N", srcstr))
        ## end if
        if not all(map(lambda v: v in atomNum, self.atom_syms)):
            raise(GErr(GErr.badatom,    # pragma: no cover
                    "Invalid atoms in list: {0}".format(self.atom_syms),
                    srcstr))
        ## end if

    ## end def __init__


    def check_geom(self, coords, atoms, tol=_DEF.GRAD_Coord_Match_Tol):
        """ Check for consistency of gradient geometry with input coords/atoms.

        The cartesian coordinates associated with a gradient object
        are considered consistent with the input `coords` if each
        component matches to within `tol`.  If
        `coords` or `atoms` vectors are passed that are of different
        length than those stored in the instance, a |False| value is
        returned, rather than an exception raised.

        Parameters
        ----------
        coords : length-3N ``np.float_``
            Vector of stacked 'lab-frame' Cartesian coordinates

        atoms  : length-N str or int
            Vector of atom symbols or atomic numbers

        tol    : float, optional
            Tolerance for acceptable deviation of each passed geometry
            coordinate  from that in the instance to still be considered
            matching. Default value is
            :data:`DEF.GRAD_Coord_Match_Tol
            <opan.const.DEF.GRAD_Coord_Match_Tol>`

        Returns
        -------
        match  : bool
            Whether input `coords` and `atoms` match those in the
            instance (|True|) or not (|False|)

        fail_type  : str
            If `match` == |False|, a string description code for the reason
            for the failed match:
            coord_dim_mismatch  : Mismatch in coordinate vector sizes
            atom_dim_mismatch   : Mismatch in atom symbol vector sizes
            coord_mismatch      : Mismatch in one or more coordinates
            atom_mismatch       : Mismatch in one or more atoms
            #DOC: Propagate info when mismatch code converted to Enum

        fail_loc   : length-3N bool or length-N bool
            np.array indicating positions of mismatch in
            either `coords` or `atoms`, depending on the value of `fail_type`.
            |True| elements indicate corresponding **MATCHING** values;
            |False| elements mark **MISMATCHES**.

        Raises
        ------
        ~exceptions.ValueError
            If ``len(coords) != 3 * len(atoms)``
        """

        # Import(s)
        from .utils import check_geom as ucg

        # Wrapper call
        result = ucg(self.geom, self.atom_syms, coords, atoms, tol=tol)

        # Return result
        return result

    ## end def check_geom

## end class SuperOpanGrad



class OrcaEngrad(SuperOpanGrad):
    """ Container for ENGRAD data generated by |orca|.

    Initialize by passing the path to the file to be loaded as the
    `path` keyword argument.

    Key information contained includes the gradient, the energy, the number
    of atoms, the geometry, and the atom IDs.  For |orca|, the precision of the
    geometry is inferior to that in an XYZ file.


    |

    **Class Variables**

    .. class:: Pat

        :func:`re.compile` patterns for data parsing.

        |

        .. attribute:: atblock

            Captures the entire block of atom ID & geometry data.

        .. attribute:: atline

            Extracts single lines from the atom ID / geometry block.

        .. attribute:: energy

            Captures the electronic energy.

        .. attribute:: gradblock

            Captures the gradient data block.

        .. attribute:: numats

            Retrieves the stand-along 'number of atoms' field.

    |

    **Instance Variables**

    atom_syms   : length-N list of str
        Uppercased atomic symbols for the atoms in the system.

    energy      : float
        Single-point energy for the geometry.

    geom        : length-3N ``np.float_``
        Vector of the atom coordinates in :math:`\\mathrm B`.

    gradient    : length-3N ``np.float_``
        Vector of the Cartesian gradient in
        :math:`\\frac{\\mathrm{E_h}}{\\mathrm B}`.

    in_str      : str
        Complete text of the ENGRAD file read in to generate the
        :class:`OrcaEngrad` instance.

    num_ats     : int
        Number of atoms in the geometry ('N')

    """

    # Imports

    # Various class-level RegEx patterns, wrapped in a class
    class Pat(object):
        # Imports
        import re as _re

        # Number of atoms -- NOT multilined. This may be a bit
        #  heavy-handed, as it will break if the .engrad file format
        #  is changed very much, but such a breakage will be a
        #  *feature*, not a bug, since it will prompt re-evaluation
        #  of the code to ensure the data is being read/parsed
        #  appropriately.
        numats = _re.compile("""
        \\#.*                   # Key text is in a comment block
        Number\\ of\\ atoms     # Find the key text
        .*\\n                   # All else from key text line, to newline
        \\#.*\\n                # Another comment line w/pound
        [ ]*                    # Series of possible spaces
        (?P<num>\\d+)           # Actual number of atoms
        .*\\n                   # All else to newline
        \\#                     # Pound to start the next line
        """, _re.I | _re.X)

        # Energy of associated geometry -- also NOT multilined.
        #  Also heavy-handed, but again this will assist in
        #  forcing review of the code if future
        #  ORCA versions change the ENGRAD file formatting
        energy = _re.compile("""
        \\#.*                               # Key text is in a comment block
        current\\ total\\ energy\\ in\\ Eh  # Key text
        .*\\n                               # Clear to newline
        \\#.*\\n                            # Blank comment line
        [ ]*                                # Series of spaces
        (?P<en>[-]?[0-9]+\\.[0-9]+)         # Energy value
        .*\\n                               # Clear to newline
        \\#                                 # Pound to start the next line
        """, _re.I | _re.X)

        # Gradient block
        gradblock = _re.compile("""
        \\#.*                           # Key text is in a comment block
        in\\ Eh/bohr                    # Key text
        .*\\n                           # Clear to newline
        \\#.*\\n                        # Blank comment line
        (?P<block>                      # Retrieving the whole grad block here
            ([ ]*[0-9.-]+[ ]*\\n)+      # Grad block entries
        )                               # Close the whole grad block
        \\#                             # Pound to signify end of block
        """, _re.I | _re.X)

        # Should probably consider going ahead and storing the indicated geometry.
        #  Would want a cross-check that yes, the geometry for which the gradient
        #  is stored DOES actually match that expected by the broader computation.
        # Should not need the same level of detailed retrieval as put into
        #  xyz, though -- just a function to confirm within some specified
        #  tolerance that a passed-in geometry is consistent with the geometry
        #  found in the ENGRAD file
        atblock = _re.compile("""
        \\#.*                           # Key text is in a comment block
        coordinates\\ in\\ Bohr         # Key text
        .*\\n                           # Clear to newline
        \\#.*                           # Blank comment line
        (?P<block>                      # Retrieving the whole geometry block
            (                           # Open group for each line
                \\n                     # Start each line with newline
                [ \\t]*                 # Leading whitespace
                ([a-z]+|\\d+)+          # Atomic symbol or number
                (                       # Group for coordinates
                    [ \\t]+             # Whitespace
                    [0-9.-]+            # Coordinate (no sci notation assumed)
                ){3}                    # Three coordinates each
                [ \\t]*                 # Optional whitespace
            )+                          # Whatever number of atoms
        )                               # Close the geometry block
        """, _re.M | _re.I | _re.X)

        # Separate Regex to retrieve the individual lines of the geometry block
        #  for subsequent parsing and storage/comparison. Not using multiline
        #  because trapping newlines seems to suffice.
        atline = _re.compile("""
        \\n                                 # Start of each line is a newline
        [ \\t]*                             # Leading whitespace
        (?P<at>([a-z]+|\\d+)+)              # Atomic symbol or number
        [ \\t]+                             # Whitespace
        (?P<c1>[0-9.-]+)                    # First coordinate
        [ \\t]+                             # Whitespace
        (?P<c2>[0-9.-]+)                    # Second coordinate
        [ \\t]+                             # Whitespace
        (?P<c3>[0-9.-]+)                    # Third coordinate
        """, _re.I | _re.X)

    ## end class Pat

    def _load(self, **kwargs):
        """ Initialize OrcaEngrad gradient object from .engrad file

        Searches indicated file for energy, geometry, gradient, and number
        of atoms and stores in the corresponding instance variables.


        Parameters
        ----------
        path : str
            Complete path to the .engrad file to be read.

        Raises
        ------
        ~opan.error.GradError
            If indicated gradient file is malformed in some fashion
        ~exceptions.IOError
            If the indicated file does not exist or cannot be read

        """

        # Imports
        from .const import CIC, atomSym, atomNum
        from .error import GradError
        from .utils import safe_cast as scast
        import numpy as np

        # Check if instantiated; complain if so
        if 'engrad_path' in dir(self):
            raise(GradError(GradError.overwrite,
                    "Cannot overwrite contents of existing OrcaEngrad", ""))
        ## end if

        # Retrieve the file target
        engrad_path = kwargs['path']

        # Open file, read contents, close stream
        with open(engrad_path,'rU') as in_fl:
            self.engrad_path = engrad_path
            self.in_str = in_fl.read()

        # Check to ensure all relevant data blocks are found
        if not self.Pat.numats.search(self.in_str):
            raise(GradError(GradError.numats,
                    "Number of atoms specification not found",
                    "ENGRAD File: {0}".format(engrad_path)))
        ## end if
        if not self.Pat.energy.search(self.in_str):
            raise(GradError(GradError.en,
                    "Energy specification not found",
                    "ENGRAD File: {0}".format(engrad_path)))
        ## end if
        if not self.Pat.gradblock.search(self.in_str):
            raise(GradError(GradError.gradblock,
                    "Gradient data block not found",
                    "ENGRAD File: {0}".format(engrad_path)))
        ## end if
        if not self.Pat.atblock.search(self.in_str):
            raise(GradError(GradError.geomblock,
                    "Geometry data block not found",
                    "ENGRAD File: {0}".format(engrad_path)))
        ## end if

        # Retrieve the number of atoms
        self.num_ats = np.int_(self.Pat.numats.search(self.in_str).group("num"))

        # Retrieve the energy
        self.energy = np.float_(self.Pat.energy.search(self.in_str).group("en"))

        # Retrieve the gradient and store numerically. Raise an error if
        #  the number of gradient elements is not equal to three times the
        #  number of atoms.
        grad_str = self.Pat.gradblock.search(self.in_str).group("block")
        if not len(grad_str.splitlines()) == 3 * self.num_ats:
            raise(GradError(GradError.gradblock,
                    "Gradient block size mismatch with number of atoms",
                    "ENGRAD File: {0}".format(engrad_path)))
        ## end if
        self.gradient = np.array(grad_str.splitlines(), dtype=np.float_)

        # Pull and store the geometry block
        geom_str = self.Pat.atblock.search(self.in_str).group("block")

        # Confirm the correct number of atoms
        if not len(self.Pat.atline.findall(geom_str)) ==  self.num_ats:
            raise(GradError(GradError.geomblock,
                    "Inconsistent number of atom coordinates in \
                    geometry block", "ENGRAD File: {0}".format(engrad_path)))
        ## end if

        # Initialize the atom symbols list
        self.atom_syms = []

        # Assemble the coordinates vector and assemble/check the element
        #  ID vector.
        # Init the atom counter and the coordinates vector
        atom_count = 0
        self.geom = np.zeros((0,), dtype=np.float_)

        # Iterate over the atom spec lines and store element and coordinates
        for line_mch in self.Pat.atline.finditer(geom_str):
            # Populate the atom symbols list; have to check for
            #  whether it's an atomic number or an element symbol
            if str.isdigit(line_mch.group("at")):
                # Atomic number
                # Check for valid number
                at_num = scast(line_mch.group("at"), np.int_)
                if not (CIC.Min_Atomic_Num <= at_num <= CIC.Max_Atomic_Num):
                    raise(GradError(GradError.geomblock,
                            "Atom #{0} is an unsupported element"
                                                        .format(atom_count),
                             "ENGRAD file: {0}".format(self.engrad_path)))
                ##end if

                # Tag on the new symbol
                self.atom_syms.append(atomSym[at_num])

            else:
                # Element symbol; store as all caps
                # Check for valid element, first by catching if the
                #  specified element string is even valid
                try:
                    at_num = atomNum[line_mch.group("at").upper()]
                except KeyError:
                    raise(GradError(GradError.geomblock,
                            "Atom #{0} is an unrecognized element"
                                                        .format(atom_count),
                            "ENGRAD file: {0}".format(engrad_path)))
                ## end try

                # Now check whether the successfully converted atomic
                #  number is in the valid range (should be a redundant check.)
                if not (CIC.Min_Atomic_Num <= at_num <= CIC.Max_Atomic_Num):
                    raise(GradError(GradError.geomblock,
                            "Atom #{0} is an unsupported element"
                                                            .format(atom_count),
                             "ENGRAD file: {0}".format(engrad_path)))
                ## end if

                # Tag on the new symbol
                self.atom_syms.append(line_mch.group("at").upper())
            ## end if

            # Append the three coordinates of the current atom to the
            #  temp coordinates vector. Unit conversion not needed since
            #  ENGRAD files report coordinates in Bohrs.
            self.geom = np.concatenate(
                    (
                        self.geom,
                        [scast(line_mch.group("c{0}".format(i)), np.float_)
                                                        for i in range(1,4)]
                     ))
        ## next line_mch

    ## end def __init__

## end def OrcaEngrad


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable")
