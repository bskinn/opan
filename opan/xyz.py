#-------------------------------------------------------------------------------
# Name:        xyz
# Purpose:     Encapsulation of parsing/handling of data from OpenBabel-type
#                molecule geometry files.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     13 Aug 2014
# Copyright:   (c) Brian Skinn 2016
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

""" Module implementing OpenBabel XYZ parsing and interpretation.

The single class :class:`~opan.xyz.OpanXYZ` imports molecular geometries
in the OpenBabel `XYZ format
<http://openbabel.org/wiki/XYZ_(format)>`__ |extlink|, with
the following variations:

 * Coordinates of any precision will be read, not just the
   :math:`10.5` specified by OpenBabel

 * Both atomic symbols and atomic numbers are valid

 * Multiple geometries/frames are supported, but the number of atoms and
   their sequence in the atoms list must be maintained in all geometries.

**Contents**

    :ref:`Class Variables <class-variables>`

    :ref:`Instance Variables <instance-variables>`

    :ref:`Methods <methods>`

        All return values from a single indicated geometry.

        :func:`~opan.xyz.OpanXYZ.geom_single` -- Entire geometry vector

        :func:`~opan.xyz.OpanXYZ.displ_single` -- Displacement
        between two atoms

        :func:`~opan.xyz.OpanXYZ.dist_single` -- Euclidean distance
        between two atoms

        :func:`~opan.xyz.OpanXYZ.angle_single` -- Spanned angle of three
        atoms (one central and two distal atoms)

        :func:`~opan.xyz.OpanXYZ.dihed_single` -- Dihedral angle
        among four atoms

    .. _toc-generators:

    :ref:`Generators <generators>`

        All yielded values are composited from an arbitrary set of geometry
        and/or atom indices; indexing with negative values is supported.

        Each parameter can either be a single index, or an iterable of
        indices. If any iterables are passed, all must be the same length R,
        and a total of R values will be returned, one
        for each set of elements across the iterables.  (This behavior
        is loosly similar to Python's :func:`zip` builtin.)  Single values are
        used as provided for all values returned.

        As an additional option, a |None| value can be passed to exactly
        one parameter, which will then be assigned the full 'natural' range
        of that parameter (``range(G)`` for `g_nums` and ``range(N)`` for
        `ats_#`).

        If the optional parameter `invalid_error` is |False| (the
        default), if :exc:`~exceptions.IndexError` or
        :exc:`~exceptions.ValueError` is raised in the course of
        calculating a given value, it is ignored and a |None| value
        is returned. If |True|, then the errors are raised as-is.

        .. note::

            For :func:`~opan.xyz.OpanXYZ.angle_iter` and
            :func:`~opan.xyz.OpanXYZ.dihed_iter`, using a |None| parameter
            with `invalid_error` == |True| is **guaranteed** to raise an
            error.

        :func:`~opan.xyz.OpanXYZ.geom_iter` -- Geometries

        :func:`~opan.xyz.OpanXYZ.displ_iter` -- Displacement vectors

        :func:`~opan.xyz.OpanXYZ.dist_iter` -- Euclidean distances

        :func:`~opan.xyz.OpanXYZ.angle_iter` -- Span angles

        :func:`~opan.xyz.OpanXYZ.dihed_iter` -- Dihedral angles


|

**Class Definition**

.. autoclass:: OpanXYZ


"""

# Imports


# Debug constant
_DEBUG = False


class OpanXYZ(object):
    """ Container for OpenBabel XYZ data.

    Initializer can be called in one of two forms::

        OpanXYZ(path='path/to/file')
        OpanXYZ(atom_syms=[{array of atoms}], coords=[{array of coordinates}])

    If `path` is specified, `atom_syms` and `coords` are ignored, and
    the instance will contain all validly
    formatted geometries present in the OpenBabel file at the indicated
    path.

    .. note::

        Certain types of improperly formatted geometry blocks, such as one with
        alphabetic characters on the 'number-of-atoms' line, may not
        raise errors during loading, but instead just result in import of fewer
        geometries/frames than expected.

    If initialized with the `atom_syms` and `coords` keyword arguments,
    the instance will contain the single geometry represented by the
    provided inputs.

    In both forms, the optional keyword argument `bohrs` can be specified, to
    indicate the units of the coordinates as Bohrs (|True|)
    or Angstroms (|False|).
    Angstrom and Bohr units are the default for the `path` and
    `atom_syms`/`coords` forms,
    respectively. The units of all coordinates stored in the instance are
    **Bohrs**.

    'N' and 'G' in the below documentation refer to the number of atoms per
    geometry and the number of geometries present in the file, respectively.
    Note that **ALL** geometries present **MUST** contain the same number
    of atoms, and the elements must all **FALL IN THE SAME SEQUENCE**
    in each geometry/frame. No error will be raised if positions of
    like atoms are swapped, but for obvious reasons this
    will almost certainly cause semantic difficulties in downstream
    computations.

    .. note::

        In |orca| '.xyz' files contain the highest
        precision geometry  information of any output (save perhaps
        the textual output generated by the program), and are stored
        in Angstrom units.


    |

    .. _class-variables:

    **Class Variables**

    .. autoattribute:: p_coords
        :annotation:

    .. autoattribute:: p_geom
        :annotation:

    .. autoattribute:: LOAD_DATA_FLAG


    |

    .. _instance-variables:

    **Instance Variables**

    Except where indicated, all |str| and |list|-of-|str| values are stored as
    :data:`LOAD_DATA_FLAG` when initialized with `atom_syms`/`coords`
    arguments.

    .. attribute:: atom_syms

        length-N |str| -- Atomic symbols for the atoms (all uppercase)

    .. attribute:: descs

        length-G |str| -- Text descriptions for each geometry
        included in a loaded file

    .. attribute:: geoms

        length-G |list| of length-3N |npfloat|_ -- Molecular
        geometry/geometries read from file or passed to `coords` argument

    .. attribute:: in_str

        |str| -- Complete contents of the input file

    .. attribute:: num_atoms

        |int| -- Number of atoms per geometry, N

    .. attribute:: num_geoms

        |int| -- Number of geometries, G

    .. attribute:: XYZ_path

        |str| -- Full path to imported OpenBabel file


    |

    .. _methods:

    **Methods**

    .. automethod:: geom_single(g_num)

    .. automethod:: displ_single(g_num, at_1, at_2)

    .. automethod:: dist_single(g_num, at_1, at_2)

    .. automethod:: angle_single(g_num, at_1, at_2, at_3)

    .. automethod:: dihed_single(g_num, at_1, at_2, at_3, at_4)


    |

    .. _generators:

    **Generators**

    .. automethod:: geom_iter(g_nums)

    .. automethod:: displ_iter(g_nums, ats_1, ats_2)

    .. automethod:: dist_iter(g_nums, ats_1, ats_2)

    .. automethod:: angle_iter(g_nums, ats_1, ats_2, ats_3)

    .. automethod:: dihed_iter(g_nums, ats_1, ats_2, ats_3, ats_4)


    """

    # Imports
    import re as _re
    from .utils.decorate import arraysqueeze as _arraysqueeze

    # Class constants
    #: |str| -- Flag for irrelevant data in `atom_syms`/`coords`
    #: initialization mode.
    LOAD_DATA_FLAG = "NOT FILE"

    # Define the RegEx to pull entire geometry blocks
    # re.M is required in order for the "^" at the start of the pattern
    #  to match more than just the first geometry block in a file
    # The "[^\\n\\w]" pattern is a weak equivalent of \\s, made necessary
    #  since \\s *includes* \\n.
    # Perhaps "[ \\t]" will suffice...
    #: |re.compile| pattern -- Retrieves full OpenBabel
    #: XYZ geometry frames/blocks
    p_geom = _re.compile("""
    ^(?P<num>\\d+)      # Integer number of atoms
    [ \\t]*\\n          # Whitespace scrub to end of line
    (?P<desc>[^\\n]+)   # Geometry description
    \\n                 # Newline to end geometry description
    (?P<coord>          # Open group for entire set of atom coords
        (                   # Open group for one coord line
            [ \\t]*?            # Any whitespace before element
            ([a-z]+|\\d+)       # Element name or atomic #
            (                   # Open group for coord value
                [ \\t]+?            # Any whitespace before coord value
                [0-9.-]+            # Decimal coordinate value
            ){3}                # Three coordinates to retrieve
            [ \\t]*\\n          # Any whitespace to end of line
        )+                  # Close one coord line group; should be >= 1
    )                   # Close group for entire set of atom coords
    """, _re.I | _re.X | _re.M)

    # Define the RegEx to pull individual lines from the coordinate block
    #  retrieved by the p_geom RegEx.
    # As w/p_geom, re.M is required for "^" to match the start of each line of
    #  the coords block.
    #: |re.compile| pattern -- Retrieves individual lines from
    #: coordinate blocks matched by :attr:`p_geom`
    p_coords = _re.compile("""
    ^[ \\t]*?               # Any whitespace before element id
    (?P<el>([a-z]+|\\d+))   # Element ID, symbol or atomic number
    [ \\t]+?                # Whitespace before first coordinate
    (?P<c1>[0-9.-]+)        # First coordinate
    [ \\t]+?                # Whitespace before second coordinate
    (?P<c2>[0-9.-]+)        # Second coordinate
    [ \\t]+?                # Whitespace before third coordinate
    (?P<c3>[0-9.-]+)        # Third coordinate
    """, _re.I | _re.X | _re.M)


    def __init__(self, **kwargs):
        """ Initialize XYZ data from file or existing geometry.

        Initializer can be called in one of two forms:

        OpanXYZ(path='path/to/file')
        OpanXYZ(atom_syms={array of atoms}, coords={array of coordinates})

        In both forms, the optional keyword argument 'bohrs' can be specified,
        to indicate the units of the coordinates as Bohrs (True) or
        Angstroms (False). Angstrom and Bohr units are the default for the
        first and second forms, respectively. The units of all coordinates
        stored in the instance are *Bohrs*, regardless of this setting.
        """

        # Two situations. One should have 'path' and, optionally, 'bohrs'
        #  specified (default here is Angstroms). The other needs an N x 1
        #  'atom_syms' vector and a 3N x 1 'coords' vector.
        if 'path' in kwargs:
            # All set for load from file. 'bohrs' defaults False here
            self._load_file(kwargs['path'], \
                    bohrs=(kwargs['bohrs'] if 'bohrs' in kwargs else False))
        else:
            # Look for the from-coordinates objects
            if all((k in kwargs) for k in ('atom_syms', 'coords')):
                # Found. Call the by-data initializer. 'bohrs' defaults to
                #  TRUE here (presume main call is from repo)
                self._load_data(kwargs['atom_syms'], \
                        kwargs['coords'], \
                        bohrs=(kwargs['bohrs'] if 'bohrs' in kwargs else True))
            else:
                # NOT found -- error!
                raise NameError("Insufficient named parameters found")
            ## end if
        ## end if

    ## end def __init__


    @_arraysqueeze(1,2)
    def _load_data(self, atom_syms, coords, bohrs=True):
        """ Internal function for making XYZ object from explicit geom data.

        Parameters
        ----------
        atom_syms
            Squeezes to array of N |str| --
            Element symbols for the XYZ. Must be valid elements as defined in
            the keys of :data:`const.atom_num <opan.const.atom_num>`.

        coords
            Squeezes to array of 3N |npfloat|_ castables --
            Coordinates for the geometry.

        bohrs
            |bool|, optional --
            Units of coordinates (default |True|)

        Raises
        ------
        ~opan.XYZError
            (typecode :attr:`~opan.error.XYZError.OVERWRITE`)
            If :class:`ORCA_XYZ` object has already been initialized.

        ~exceptions.ValueError
            If atom_syms & coords dimensions are incompatible

        ~exceptions.ValueError
            If type of `atom_syms` and/or `coords` is invalid

        """

        # Imports
        import numpy as np
        from .const import atom_num
        from .error import XYZError

        # Gripe if already initialized
        if 'geoms' in dir(self):
            raise XYZError(XYZError.OVERWRITE,
                    "Cannot overwrite contents of existing OpanXYZ", "")
        ## end if

        # Check and store dimensions
        if not len(coords.shape) == 1:
            raise ValueError("Coordinates are not a vector")
        ## end if
        if not len(atom_syms.shape) == 1:
            raise ValueError("Atom symbols are not a simple list")
        ## end if
        if not coords.shape[0] == 3 * atom_syms.shape[0]:
            raise ValueError("len(coords) != 3 * len(atom_syms)")
        ## end if

        # Proof the atoms list
        if not all( (atom_syms[i].upper() in atom_num)
                                for i in range(atom_syms.shape[0]) ):
            # Invalid atoms specified
            raise ValueError("Invalid atoms specified: {0}".format(
                    [(j, atom_syms[j]) for j in
                        (i for (i, valid) in
                            enumerate(map(lambda k: k in atom_num, atom_syms))
                            if not valid
                        )
                    ] ))
        ## end if

        # Ensure the geometry is all numeric
        if not all(map(np.isreal, coords)):
            raise ValueError("All coordinates must be real numeric")
        ## end if

        # Store the number of atoms. Only one geometry. Standard string
        #  content for things only relevant to file load.
        self.num_atoms = atom_syms.shape[0]
        self.num_geoms = 1
        self.in_str = self.LOAD_DATA_FLAG
        self.descs = np.array([self.LOAD_DATA_FLAG])
        self.XYZ_path = self.LOAD_DATA_FLAG

        # Store the atoms as vector
        self.atom_syms = list(map(str.upper, list(atom_syms)))

        # Store the single geometry by bracketing with an array
        self.geoms = [coords / (1.0 if bohrs else PHYS.ANG_PER_BOHR)]

    ## end def _load_data


    def _load_file(self, XYZ_path, bohrs=False):
        """ Initialize OpanXYZ geometry object from OpenBabel file

        Import of an arbitrary number of multiple geometries from an OpenBabel
        file. All geometries must have the same number and type of atoms,
        and the ordering of atom types must be retained throughout.

        Parameters
        ----------
        XYZ_path
            |str| --
            Complete path to the .xyz / .trj file to be read

        bohrs
            |bool|, optional --
            Flag for whether coordinates are interpreted as in Bohrs or
            Angstroms. Default is |False| (Angstroms).

        Raises
        ------
        ~opan.error.XYZError
            (typecode :attr:`~opan.error.XYZError.XYZFILE`)
            If indicated geometry file is malformed in some fashion

        ~exceptions.IOError
            If the indicated file does not exist or cannot be read

        """

        # Imports
        import numpy as np
        from .const import CIC, PHYS, atom_num, atom_sym
        from .error import XYZError
        from .utils import safe_cast as scast

        # Complain if already initialized
        if 'geoms' in dir(self):
            raise XYZError(XYZError.OVERWRITE,
                    "Cannot overwrite contents of existing OpanXYZ", "")
        ## end if

        # Open file, read contents, close stream
        # No particular exception handling; that will be the responsibility
        #  of the calling routine.
        # May want to add a check for text file format; but, really, a check
        #  for no-found-geometry-information will likely cover such a case.
        with open(XYZ_path,'rU') as in_fl:
            self.XYZ_path = XYZ_path #: Path of file
            self.in_str = in_fl.read()

        # Check to ensure at least one geom match; else raise some sort of
        #  error
        if not OpanXYZ.p_geom.search(self.in_str):
            raise XYZError(XYZError.XYZFILE,
                    "No valid geometry found",
                    "XYZ file: " + XYZ_path)
        ## end if

        # Store the number of atoms. XYZ files with multiple geometries will
        #  be required to have the same number of atoms in all geometries.
        #  This is fine for OpenAnharmonic, but may be highly undesirable
        #  in other applications.  Beware.
        # If the .match() call doesn't have a group attribute, catch the
        #  resulting AttributeError and throw an XYZError. Can't think of any
        #  other situation that would throw an AttributeError here, so probably
        #  I don't have to check for a subtype of AttributeError...
        # Can't think of a need for an else or a finally here
        try:
            self.num_atoms = scast(OpanXYZ.p_geom.match(self.in_str)
                                            .group("num"), np.int_)
        except AttributeError:
            raise XYZError(XYZError.XYZFILE,
                    "No geometry block found at start of file",
                    "XYZ file: " + XYZ_path)
        ## end try

        # Initialize the description vector and geometry array
        # CANNOT use matrix type for the geometry since it's 3D.
        # Retriever function for geometry should convert the given slice
        #  to a vector as an np.matrix, though.
        # Geometry is constructed as a length-3N vector of
        #  coordinates, where atom 1's x/y/z coordinates are sequential,
        #  then atom 2's x/y/z coordinates are sequential, etc.
        # Multiple geometries in a single XYZ file are stacked in a simple
        #  array.
        self.descs = []
        self.geoms = []

        # Initialize the atom symbols vector
        self.atom_syms = []

        # Define a counter for the number of geometries. Store as a persistent
        #  instance variable because it will be useful later.
        self.num_geoms = 0

        # Loop over the geometry blocks found in the input file
        for mch in OpanXYZ.p_geom.finditer(self.in_str):
            # Check that the number of atoms is consistent with the spec
            #  found in the first geometry block
            if not scast(mch.group("num"), np.int_) == self.num_atoms:
                raise XYZError(XYZError.XYZFILE,
                        "Non-constant number of atoms in multiple geometry",
                        "XYZ file: " + XYZ_path)
            ## end if

            # Store the description for the current geometry
            self.descs.append(mch.group("desc"))

            # Assemble the coordinates vector and assemble/check the element
            #  ID vector.
            # Reset the atom counter and the coordinates vector
            atom_count = 0
            coord_vec = np.empty((0,),dtype=np.float_)
            for line_mch in OpanXYZ.p_coords.finditer(mch.group("coord")):
                # Check for whether element list has been fully populated
                if len(self.atom_syms) < self.num_atoms:
                    # If not, continue populating it; have to check for
                    #  whether it's an atomic number or an element symbol
                    if line_mch.group("el").isdigit():
                        # Atomic number
                        # Check for valid number
                        at_num = scast(line_mch.group("el"), np.int_)
                        if not (CIC.MIN_ATOMIC_NUM <= at_num
                                                    <= CIC.MAX_ATOMIC_NUM):
                            raise XYZError(XYZError.XYZFILE,
                                    "Geometry #{0}, atom #{1} is an \
                                        unsupported element"
                                        .format(self.num_geoms, atom_count),
                                     "XYZ file: {0}".format(XYZ_path))
                        ##end if

                        # Tag on the new symbol
                        self.atom_syms.append(atom_sym[at_num])
                    else:
                        # Element symbol; store as all caps
                        # Check for valid element, first by catching if the
                        #  specified element string is even valid
                        try:
                            at_num = atom_num[line_mch.group("el").upper()]
                        except KeyError:
                            raise XYZError(XYZError.XYZFILE,
                                    "Geometry #{0}, atom #{1} is an \
                                        unrecognized element"
                                        .format(self.num_geoms, atom_count),
                                     "XYZ file: {0}".format(XYZ_path))
                        ## end try

                        # Tag on the new symbol
                        self.atom_syms.append(line_mch.group("el").upper())
                    ## end if
                else: # List is fully populated
                    # If so, confirm that the element specified at the current
                    #  line matches that of the first geometry of the file
                    # Have to check whether it's an atomic number or symbol.
                    if line_mch.group("el").isdigit():
                        # Atomic number; don't need safe cast; must trap for
                        #  a bad atomic number
                        at_num = scast(line_mch.group("el"), np.int_)
                        if not (CIC.MIN_ATOMIC_NUM <= at_num
                                                    <= CIC.MAX_ATOMIC_NUM):
                            raise XYZError(XYZError.XYZFILE,
                                    "Geometry #{0}, atom #{1} is an \
                                        unsupported element"
                                        .format(self.num_geoms, atom_count),
                                     "XYZ file: {0}".format(XYZ_path))
                        ## end if
                        if not atom_sym[at_num] == self.atom_syms[atom_count]:
                            raise XYZError(XYZError.XYZFILE,
                                    "Geometry #{0}, atom #{1} is inconsistent \
                                        with geometry #0"
                                        .format(self.num_geoms, atom_count),
                                     "XYZ file: {0}".format(XYZ_path))
                        ## end if

                    else:
                        # Element symbol
                        # Check for valid element, first by catching if the
                        #  specified element string is even valid
                        try:
                            at_num = atom_num[line_mch.group("el").upper()]
                        except KeyError:
                            raise XYZError(XYZError.XYZFILE,
                                    "Geometry #{0}, atom #{1} is an \
                                        unrecognized element"
                                        .format(self.num_geoms, atom_count),
                                     "XYZ file: {0}".format(XYZ_path))
                        ## end try
                        # Confirm symbol matches the initial geometry
                        if not line_mch.group("el").upper() == \
                                                self.atom_syms[atom_count]:
                            raise XYZError(XYZError.XYZFILE,
                                    "Geometry #{0}, atom #{1} is inconsistent \
                                        with geometry #0"
                                        .format(self.num_geoms, atom_count),
                                    "XYZ file: {0}".format(XYZ_path))
                        ## end if
                    ## end if
                ## end if

                # Append the three coordinates of the current atom to the
                #  temp coordinates vector, converting to Bohrs if indicated.
                # Working in Bohrs is desired because they are atomic units
                #  and thus are part of the internal unit definition of the
                #  Hartree.
                for c_str in range(1,4):
                    coord_vec = np.concatenate(
                            (coord_vec, [
                        scast(line_mch.group("c{0}".format(c_str)), np.float_) /
                        (1.0 if bohrs else PHYS.ANG_PER_BOHR)
                                        ]), axis=0)
                ## next c_str

                # Increment the atom_count for the atom number
                atom_count += 1

            ## next line_mch

            # Confirm that number of imported coordinates matches the
            #  number expected from self.num_atoms
            if not coord_vec.shape[0] == 3*self.num_atoms:
                raise XYZError(XYZError.XYZFILE,
                        "Geometry #{0} atom count is inconsistent"
                            .format(self.num_geoms),
                        "XYZ file: {0}".format(XYZ_path))
            ## end if

            # Assemble the coordinates vector into the actual coordinates
            #  stack for the XYZ
            self.geoms.append(coord_vec)

            # Increment the count of the number of geometries. Once the
            #  'mch' iteration is completed, this will accurately reflect
            #  the number of geometries read from the file.
            self.num_geoms += 1

        ## next mch

    ## end def _load_file


    def geom_single(self, g_num):
        """ Retrieve a single geometry.

        The atom coordinates are returned with each atom's
        x/y/z coordinates grouped together::

            [A1x, A1y, A1z, A2x, A2y, A2z, ...]

        An alternate method to achieve the same effect as this
        function is by simply
        indexing into :attr:`~opan.xyz.OpanXYZ.geoms`::

            >>> x = opan.xyz.OpanXYZ(path='...')
            >>> x.geom_single(g_num)    # One way to do it
            >>> x.geoms[g_num]          # Another way

        Parameters
        ----------
        g_num
            |int| --
            Index of the desired geometry

        Returns
        -------
        geom
            length-3N |npfloat|_ --
            Vector of the atomic coordinates for the geometry indicated
            by `g_num`

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` is provided
        """

        # Just return the appropriate geometry vector
        geom = self.geoms[g_num]
        return geom

    ## end def geom_single


    def geom_iter(self, g_nums):
        """Iterator over a subset of geometries.

        The indices of the geometries to be returned are indicated by an
        iterable of |int|\\ s passed as `g_nums`.

        As with :meth:`geom_single`, each geometry is returned as a
        length-3N |npfloat|_ with each atom's x/y/z coordinates
        grouped together::

            [A1x, A1y, A1z, A2x, A2y, A2z, ...]

        In order to use NumPy `slicing or advanced indexing
        <http://docs.scipy.org/doc/numpy-1.10.0/reference/
        arrays.indexing.html>`__, :data:`geoms` must first be
        explicitly converted to |nparray|, e.g.::

            >>> x = opan.xyz.OpanXYZ(path='...')
            >>> np.array(x.geoms)[[2,6,9]]

        Parameters
        ----------
        g_nums
            length-R iterable of |int| --
            Indices of the desired geometries

        Yields
        ------
        geom
            length-3N |npfloat|_ --
            Vectors of the atomic coordinates for each geometry
            indicated in `g_nums`

        Raises
        ------
        ~exceptions.IndexError
            If an item in `g_nums` is invalid (out of range)

        """
        # Using the custom coded pack_tups to not have to care whether the
        #  input is iterable
        from .utils import pack_tups

        vals = pack_tups(g_nums)
        for val in vals:
            yield self.geom_single(val[0])

    ## end def geom_iter


    def dist_single(self, g_num, at_1, at_2):
        """ Distance between two atoms.

        Parameters
        ----------
        g_num
            |int| -- Index of the desired geometry

        at_1
            |int| -- Index of the first atom

        at_2
            |int| -- Index of the second atom

        Returns
        -------
        dist
            |npfloat|_ --
            Distance in Bohrs between `at_1` and `at_2` from
            geometry `g_num`

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided

        """

        # Import used math library function(s)
        import numpy as np
        from scipy import linalg as spla
        from .utils import safe_cast as scast

        # The below errors are explicitly thrown since values are multiplied by
        #  three when they are used as an index and thus give non-intuitive
        #  errors in subsequent code.
        # Complain if at_1 is invalid
        if not (-self.num_atoms <= at_1 < self.num_atoms):
            raise IndexError("Invalid index for 'at_1' ({0})".format(at_1))

        # Complain if at_2 is invalid
        if not (-self.num_atoms <= at_2 < self.num_atoms):
            raise IndexError("Invalid index for 'at_2' ({0})".format(at_2))

        # Should never be necessary (save for badly erroneous calling code),
        #  but coerce at_1 and at_2 to their floor() values.  This is again
        #  needed since they are multiplied by three in the index expresssions
        #  below, and can cause funny behavior when truncated by the indexing
        at_1 = scast(np.floor(at_1), np.int_)
        at_2 = scast(np.floor(at_2), np.int_)

        # Calculate the interatomic distance and return. Return identically
        #  zero if the indices are equal
        if at_1 == at_2:
            dist = 0.0
        else:
            dist = scast(
                    spla.norm(self.displ_single(g_num, at_1, at_2)),
                            np.float_)
        ## end if

        return dist

    ## end def dist_single


    def dist_iter(self, g_nums, ats_1, ats_2, invalid_error=False):
        """ Iterator over selected interatomic distances.

        Distances are in Bohrs as with :meth:`dist_single`.

        See `above <toc-generators_>`_ for more information on
        calling options.

        Parameters
        ----------
        g_nums
            |int| or length-R iterable |int| or |None| --
            Index/indices of the desired geometry/geometries

        ats_1
            |int| or iterable |int| or |None| --
            Index/indices of the first atom(s)

        ats_2
            |int| or iterable |int| or |None| --
            Index/indices of the second atom(s)

        invalid_error
            |bool|, optional --
            If |False| (the default), |None| values are returned for
            results corresponding to invalid indices. If |True|,
            exceptions are raised per normal.

        Yields
        ------
        dist
            |npfloat|_ --
            Interatomic distance in Bohrs between each atom pair of
            `ats_1` and `ats_2` from the corresponding geometries
            of `g_nums`.

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided.

        ~exceptions.ValueError
            If all iterable objects are not the same length.

        """

        # Imports
        import numpy as np
        from .utils import pack_tups

        # Print the function inputs if debug mode is on
        if _DEBUG:  # pragma: no cover
            print("g_nums = {0}".format(g_nums))
            print("ats_1 = {0}".format(ats_1))
            print("ats_2 = {0}".format(ats_2))
        ## end if

        # Perform the None substitution
        arglist = self._none_subst(g_nums, ats_1, ats_2)

        # Expand/pack the tuples from the inputs
        tups = pack_tups(*arglist)

        # Dump the results if debug mode is on
        if _DEBUG:  # pragma: no cover
            print(tups)
        ## end if

        # Construct the generator using the packed tuples. If 'None' expansion
        #  was used, return None for any invalid indices instead of raising
        #  an exception.
        for tup in tups:
            yield self._iter_return(tup, self.dist_single, invalid_error)
        ## next tup

    ## end def dist_iter


    def angle_single(self, g_num, at_1, at_2, at_3):
        """ Spanning angle among three atoms.

        The indices `at_1` and `at_3` can be the same (yielding a
        trivial zero angle), but `at_2` must be different from
        both `at_1` and `at_3`.

        Parameters
        ----------
        g_num
            |int| --
            Index of the desired geometry

        at_1
            |int| --
            Index of the first atom

        at_2
            |int| --
            Index of the second atom

        at_3
            |int| --
            Index of the third atom

        Returns
        -------
        angle
            |npfloat|_ --
            Spanning angle in degrees between `at_1`-`at_2`-`at_3`, from
            geometry `g_num`

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided

        ~exceptions.ValueError
            If `at_2` is equal to either `at_1` or `at_3`

        """

        # Imports
        import numpy as np
        from .utils import safe_cast as scast
        from .utils.vector import vec_angle

        # The below errors are explicitly thrown since they are multiplied by
        #  three when they are used as an index and thus give non-intuitive
        #  errors in later code.
        # Complain if at_1 is invalid
        if not(-self.num_atoms <= at_1 < self.num_atoms):
            raise IndexError("Invalid index for 'at_1' ({0})".format(at_1))

        # Complain if at_2 is invalid
        if not(-self.num_atoms <= at_2 < self.num_atoms):
            raise IndexError("Invalid index for 'at_2' ({0})".format(at_2))

        # Complain if at_3 is invalid
        if not(-self.num_atoms <= at_3 < self.num_atoms):
            raise IndexError("Invalid index for 'at_3' ({0})".format(at_3))

        # Should never be necessary (save for badly erroneous calling code),
        #  but coerce the at_x to their floor() values.  This is again
        #  needed since they are multiplied by three in the index expresssions
        #  below, and can cause funny behavior when truncated by the indexing
        at_1 = scast(np.floor(at_1), np.int_)
        at_2 = scast(np.floor(at_2), np.int_)
        at_3 = scast(np.floor(at_3), np.int_)

        # Complain if at_2 is equal to either at_1 or at_3.  Must factor in
        #  the possibility of negative indexing via modulo arithmetic.
        if (at_2 % self.num_atoms) == (at_1 % self.num_atoms):
            raise ValueError("'at_1' and 'at_2' must be different")
        if (at_2 % self.num_atoms) == (at_3 % self.num_atoms):
            raise ValueError("'at_2' and 'at_3' must be different")

        # Trivial return if at_1 and at_3 are the same
        if (at_1 % self.num_atoms) == (at_3 % self.num_atoms):
            # Angle is identically zero in this case
            return 0.0
        ## end if

        # Store the displacement vectors from at_2 to at_1 and to at_3
        # The np.float64 type should be retained through the displ_single call.
        vec_2_1 = self.displ_single(g_num, at_2, at_1)
        vec_2_3 = self.displ_single(g_num, at_2, at_3)

        # Compute and return the calculated angle, in degrees
        # v1 {dot} v2 == |v1||v2| * cos(theta)
        angle = vec_angle(vec_2_1, vec_2_3)
        return angle

    ## end def angle_single


    def angle_iter(self, g_nums, ats_1, ats_2, ats_3, invalid_error=False):
        """ Iterator over selected atomic angles.

        Angles are in degrees as with :meth:`angle_single`.

        See `above <toc-generators_>`_ for more information on
        calling options.

        Parameters
        ----------
        g_nums
            |int| or iterable |int| or |None| --
            Index of the desired geometry

        ats_1
            |int| or iterable |int| or |None| --
            Index of the first atom

        ats_2
            |int| or iterable |int| or |None| --
            Index of the second atom

        ats_3
            |int| or iterable |int| or |None| --
            Index of the third atom

        invalid_error
            |bool|, optional --
            If |False| (the default), |None| values are returned for
            results corresponding to invalid indices. If |True|,
            exceptions are raised per normal.

        Yields
        ------
        angle
            |npfloat|_ --
            Spanning angles in degrees between corresponding |br|
            `ats_1`-`ats_2`-`ats_3`, from geometry/geometries `g_nums`

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided.

        ~exceptions.ValueError
            If all iterable objects are not the same length.

        ~exceptions.ValueError
            If any `ats_2` element is equal to either the corresponding `ats_1`
            or `ats_3` element.

        """
        # Suitability of ats_n indices will be checked within the
        #  self.angle_single() calls and thus no check is needed here.

        # Import the tuple-generating function
        from .utils import pack_tups

        # Print the function inputs if debug mode is on
        if _DEBUG:   # pragma: no cover
            print("g_nums = {0}".format(g_nums))
            print("ats_1 = {0}".format(ats_1))
            print("ats_2 = {0}".format(ats_2))
            print("ats_3 = {0}".format(ats_3))
        ## end if

        # Perform the None substitution
        arglist = self._none_subst(g_nums, ats_1, ats_2, ats_3)

        # Expand/pack the tuples from the inputs
        tups = pack_tups(*arglist)

        # Dump the results if debug mode is on
        if _DEBUG:  # pragma: no cover
            print(tups)
        ## end if

        # Construct the generator using the packed tuples.
        for tup in tups:
            if _DEBUG: # pragma: no cover
                print(tup)
            ## end if

            yield self._iter_return(tup, self.angle_single, invalid_error)
        ## next tup

    ## end def angle_iter


    def dihed_single(self, g_num, at_1, at_2, at_3, at_4):
        """ Dihedral/out-of-plane angle among four atoms.

        Returns the out-of-plane angle among four atoms from geometry
        `g_num`, in degrees.  The reference plane
        is spanned by `at_1`, `at_2` and `at_3`. The out-of-plane angle is
        defined such that a positive angle represents a counter-clockwise
        rotation of the projected `at_3`\\ :math:`\\rightarrow`\\ `at_4`
        vector with respect to the
        reference plane when looking from `at_3` toward `at_2`.
        Zero rotation corresponds to occlusion of `at_1` and `at_4`;
        that is, the case where
        the respective rejections of `at_1`
        :math:`\\rightarrow`\\ `at_2` and
        `at_3`\\ :math:`\\rightarrow`\\ `at_4` onto
        `at_2`\\ :math:`\\rightarrow`\\ `at_3`
        are ANTI-PARALLEL.

        *#DOC: Pull the above to User Guide eventually, with figures.*

        All four atom indices must be distinct. Both of the atom trios 1-2-3
        and 2-3-4 must be sufficiently nonlinear, as diagnosed by a bend
        angle different from 0 or 180 degrees by at least
        :data:`PRM.NON_PARALLEL_TOL <opan.const.PRM.NON_PARALLEL_TOL>`.

        Parameters
        ----------
        g_num
            |int| -- Index of the desired geometry

        at_1
            |int| -- Index of the first atom

        at_2
            |int| -- Index of the second atom

        at_3
            |int| -- Index of the third atom

        at_4
            |int| -- Index of the fourth atom

        Returns
        -------
        dihed
            |npfloat|_ --
            Out-of-plane/dihedral angle in degrees for the indicated `at_#`,
            drawn from geometry `g_num`

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided

        ~exceptions.ValueError
            If any indices `at_#` are equal

        ~opan.error.XYZError
            (typecode :data:`~opan.error.XYZError.DIHED`) If either
            of the atom trios (1-2-3 or 2-3-4) is too close to linearity

        """
        # library imports
        import numpy as np
        from scipy import linalg as spla
        from .utils.vector import ortho_basis, rej, vec_angle
        from .utils import safe_cast as scast
        from .error import XYZError
        from .const import PRM

        # The below errors are explicitly checked and raised since the indices
        #  are multiplied by three when they are used as an index
        #  and thus give non-intuitive errors in later code.
        # Complain if at_1 is invalid
        if not(-self.num_atoms <= at_1 < self.num_atoms):
            raise IndexError("Invalid index for 'at_1' ({0})".format(at_1))

        # Complain if at_2 is invalid
        if not(-self.num_atoms <= at_2 < self.num_atoms):
            raise IndexError("Invalid index for 'at_2' ({0})".format(at_2))

        # Complain if at_3 is invalid
        if not(-self.num_atoms <= at_3 < self.num_atoms):
            raise IndexError("Invalid index for 'at_3' ({0})".format(at_3))

        # Complain if at_4 is invalid
        if not(-self.num_atoms <= at_4 < self.num_atoms):
            raise IndexError("Invalid index for 'at_4' ({0})".format(at_4))

        # Should never be necessary (save for badly erroneous calling code),
        #  but coerce the at_x to their floor() values.  This is again
        #  needed since they are multiplied by three in the index expresssions
        #  below, and can cause funny behavior when truncated by the indexing
        at_1 = scast(np.floor(at_1), np.int_)
        at_2 = scast(np.floor(at_2), np.int_)
        at_3 = scast(np.floor(at_3), np.int_)
        at_4 = scast(np.floor(at_4), np.int_)

        # Proofread the atom numbers. Performed by double-iterative scan of
        #  the atom numbers, converting the index equality test results to
        #  ints and summing the results.  Since each ats_n is not compared to
        #  itself, a sum of zero should diagnose the required mutually
        #  nonidentical indices.
        #
        # Pile the atom indices into a vector
        ats = [at_1, at_2, at_3, at_4]

        # Scan over the vector of indices pairwise without repetition, and
        #  without examining for at_i == at_i (which is trivially and always
        #  True).  Store the outcomes as integers (True == 1; False == 0)
        ats_test = [int(ats[x] == ats[y]) for x in range(4) \
                                        for y in range(x+1,4)]

        # For a proper set of indices, the sum over ats_test will be zero.
        if sum(ats_test) > 0:
            # Improper set of indices; at least one pair is duplicated.
            #  Collate the duplicative pairings and raise ValueError.
            #  np.triu_indices generates index pairs in the same sequence as
            #  the above double iteration over ats, but as a list of two
            #  np.arrays.  column_stack puts them together as column vectors,
            #  allowing the conditional iteration over x to select only those
            #  index pairs that correspond to duplicated indices.  The
            #  resulting filtered pairs are converted to tuples for concise
            #  formatting in the output.
            ats_pairs = [tuple(np.column_stack(np.triu_indices(4,1))[x])
                                        for x in range(6) if ats_test[x] == 1]
            raise ValueError("Duplicate atom indices: {0}".format(ats_pairs))
        ## end if

        # Check to ensure non-collinearity of the 1-2-3 and 2-3-4 atom trios
        for idx in range(2):
            # Store the relevant angle
            ang = self.angle_single(g_num, [at_2, at_3][idx],
                                         [at_1, at_2][idx],
                                         [at_3, at_4][idx])

            # Check for whether angle is too close to zero or 180 degrees
            if np.min([ang, 180.0 - ang]) < PRM.NON_PARALLEL_TOL:
                # Too close; raise error
                raise XYZError(XYZError.DIHED,
                        "Angle {0} is insufficiently nonlinear"
                            .format([(at_2, at_1, at_3),
                            (at_3, at_2, at_4)][idx]),
                        "XYZ file: {0}".format(self.XYZ_path))
            ## end if
        ## next idx

        # Store normalized atomic displacement vector at_2 --> at_3 as that
        #  defining the projection plane
        plane_norm = self.displ_single(g_num, at_2, at_3)
        plane_norm /= spla.norm(plane_norm)

        # Retrieve the orthonormal basis in the projection plane, with the
        #  first vector being the normalized projection of the at_1 --> at_2
        #  displacement onto that plane
        on1, on2 = ortho_basis(plane_norm, \
                            self.displ_single(g_num, at_1, at_2))

        # Project the at_3 --> at_4 displacement onto the plane
        #
        # Retrieve the "back-side" displacement vector
        back_vec = self.displ_single(g_num, at_3, at_4)

        # Project onto the plane by subtracting out the plane_norm projection
        #  and re-normalize
        back_vec = rej(back_vec, plane_norm)
        back_vec /= spla.norm(back_vec)

        # Calculate the absolute value of the departure of the dihedral/
        #  out-of-plane angle from 180 degrees as derived from the dot-product
        #  of on1 and back_vec. Both should be normalized at this point, so
        #  the calculation is straightforward
        dihed = vec_angle(back_vec, on1)

        # Given the handedness of the spanning vectors provided by ortho_basis,
        #  the sign of the dihed departure is that of the dot product
        #  of back_vec and on2.
        dihed *= np.sign(np.dot(back_vec, on2))

        # Conversion to the stated typical definition of a dihedral now
        #  requires addition of 180 degrees.
        dihed += 180.0

        # Should be set to return the value
        return dihed

    ## end dihed_single


    def dihed_iter(self, g_nums, ats_1, ats_2, ats_3, ats_4, \
                                                    invalid_error=False):
        """ Iterator over selected dihedral angles.

        Angles are in degrees as with :meth:`dihed_single`.

        See `above <toc-generators_>`_ for more information on
        calling options.


        Parameters
        ----------
        g_nums
            |int| or iterable |int| or |None| --
            Indices of the desired geometry

        ats_1
            |int| or iterable |int| or |None| --
            Indices of the first atoms

        ats_2
            |int| or iterable |int| or |None| --
            Indices of the second atoms

        ats_3
            |int| or iterable |int| or |None| --
            Indices of the third atoms

        ats_4
            |int| or iterable |int| or |None| --
            Indices of the fourth atoms

        invalid_error
            |bool|, optional --
            If |False| (the default), |None| values are returned for
            results corresponding to invalid indices. If |True|,
            exceptions are raised per normal.

        Yields
        ------
        dihed
            |npfloat|_ --
            Out-of-plane/dihedral angles in degrees for the indicated
            atom sets `ats_1`-`ats_2`-`ats_3`-`ats_4`, drawn from
            the respective `g_nums`.

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided.

        ~exceptions.ValueError
            If all iterable objects are not the same length.

        ~exceptions.ValueError
            If any corresponding `ats_#` indices are equal.

        ~opan.error.XYZError
            (typecode :data:`~opan.error.XYZError.DIHED`) If either
            of the atom trios (1-2-3 or
            2-3-4) is too close to linearity for any group of `ats_#`

        """
        # Suitability of ats_n indices will be checked within the
        #  self.dihed_single() calls and thus no check is needed here.

        # Import the tuple-generating function
        from .utils import pack_tups

        # Print the function inputs if debug mode is on
        if _DEBUG:   # pragma: no cover
            print("g_nums = {0}".format(g_nums))
            print("ats_1 = {0}".format(ats_1))
            print("ats_2 = {0}".format(ats_2))
            print("ats_3 = {0}".format(ats_3))
            print("ats_4 = {0}".format(ats_4))
        ## end if

        # Perform the None substitution
        arglist = self._none_subst(g_nums, ats_1, ats_2, ats_3, ats_4)

        # Expand/pack the tuples from the inputs
        tups = pack_tups(*arglist)

        # Dump the results if debug mode is on
        if _DEBUG:   # pragma: no cover
            print(tups)
        ## end if

        # Construct the generator using the packed tuples.
        for tup in tups:
            yield self._iter_return(tup, self.dihed_single, invalid_error)
        ## next tup

    ## end def dihed_iter


    def displ_single(self, g_num, at_1, at_2):
        """ Displacement vector between two atoms.

        Returns the displacement vector pointing from `at_1`
        toward `at_2` from geometry `g_num`.
        If `at_1` == `at_2` a strict zero vector is returned.

        Displacement vector is returned in units of Bohrs.

        Parameters
        ----------
        g_num
            |int| -- Index of the desired geometry

        at_1
            |int| -- Index of the first atom

        at_2
            |int| -- Index of the second atom

        Returns
        -------
        displ
            length-3 |npfloat|_ --
            Displacement vector from `at_1` to `at_2`

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided

        """

        # Library imports
        import numpy as np
        from .utils import safe_cast as scast

        # The below errors are explicitly thrown since they are multiplied by
        #  three when they are used as an index and thus give non-intuitive
        #  errors.
        # Complain if at_1 is invalid
        if not (-self.num_atoms <= at_1 < self.num_atoms):
            raise IndexError("Invalid index for 'at_1' ({0})".format(at_1))

        # Complain if at_2 is invalid
        if not (-self.num_atoms <= at_2 < self.num_atoms):
            raise IndexError("Invalid index for 'at_2' ({0})".format(at_2))

        # Should never be necessary (save for badly erroneous calling code),
        #  but coerce at_1 and at_2 to their floor() values.  This is again
        #  needed since they are multiplied by three in the index expresssions
        #  below, and can cause funny behavior when truncated by the indexing
        at_1 = scast(np.floor(at_1), np.int_)
        at_2 = scast(np.floor(at_2), np.int_)

        # If the atom indices are the same, return trivial zero vector
        if (at_1 % self.num_atoms) == (at_2 % self.num_atoms):
            return np.array([0.0, 0.0, 0.0])
        ## end if

        # Retrieve the geometry; np.float_ type should be retained
        g = self.geom_single(g_num)

        # Calculate the displacement vector and return
        displ = np.array([ g[i + 3*at_2] - g[i + 3*at_1] for i in range(3) ])

        # Return the displacement vector
        return displ

    ## end def displ_single


    def displ_iter(self, g_nums, ats_1, ats_2, invalid_error=False):
        """ Iterator over indicated displacement vectors.

        Displacements are in Bohrs as with :meth:`displ_single`.

        See `above <toc-generators_>`_ for more information on
        calling options.

        Parameters
        ----------
        g_nums
            |int| or length-R iterable |int| or |None| --
            Index/indices of the desired geometry/geometries

        ats_1
            |int| or length-R iterable |int| or |None| --
            Index/indices of the first atom(s)

        ats_2
            |int| or length-R iterable |int| or |None| --
            Index/indices of the second atom(s)

        invalid_error
            |bool|, optional --
            If |False| (the default), |None| values are returned for
            results corresponding to invalid indices. If |True|,
            exceptions are raised per normal.

        Yields
        ------
        displ
            |npfloat|_ --
            Displacement vector in Bohrs between each atom pair of |br|
            `ats_1` :math:`\\rightarrow` `ats_2` from the corresponding
            geometries of `g_nums`.

        Raises
        ------
        ~exceptions.IndexError
            If an invalid (out-of-range) `g_num` or `at_#` is provided.

        ~exceptions.ValueError
            If all iterable objects are not the same length.

        """

        # Import the tuple-generating function
        from .utils import pack_tups

        # Print the function inputs if debug mode is on
        if _DEBUG:  # pragma: no cover
            print("g_nums = {0}".format(g_nums))
            print("ats_1 = {0}".format(ats_1))
            print("ats_2 = {0}".format(ats_2))
        ## end if

        # Perform the None substitution
        arglist = self._none_subst(g_nums, ats_1, ats_2)

        # Expand/pack the tuples from the inputs
        tups = pack_tups(*arglist)

        # Dump the results if debug mode is on
        if _DEBUG:  # pragma: no cover
            print(tups)
        ## end if

        # Construct the generator using the packed tuples.
        for tup in tups:
            yield self._iter_return(tup, self.displ_single, invalid_error)
        ## next tup

    ## end def displ_iter


    def _none_subst(self, *args):
        """ Helper function to insert full ranges for |None| for X_iter methods.

        Custom method, specifically tailored, taking in the arguments from
        an X_iter method and performing the replacement of |None| after
        error-checking the arguments for a max of one |None| value, and ensuring
        that if a |None| is present, no other non-|str| iterables are present.

        Parameters
        ----------
        args : 3-5 arguments of |int| or iterable |int|, or |None|
            First argument is always the indices for the geometries; all
            following are for the atoms in sequence as required for the
            particular :samp:`{x}_iter` method

        Returns
        -------
        arglist     : 3-5 arguments, matching input params
            Argument list, with |None| substituted if validly present

        Raises
        ------
        ~exceptions.ValueError  : If more than one |None| argument is present

        ~exceptions.ValueError  : If an arg is non-|str| iterable when one
        |None| is present
        """

        # Imports
        import numpy as np

        # Initialize argument list return value, and as None not found
        arglist = [a for a in args]
        none_found = False

        # Check for None values
        none_vals = list(map(lambda e: isinstance(e, type(None)), arglist))

        # Error if more than one None; handle if exactly one; pass through if
        #  none.
        if np.count_nonzero(none_vals) > 1:
            raise ValueError("Multiple 'None' values [indices {0}] not \
                    supported".format(tuple(np.nonzero(none_vals)[0])))
        elif np.count_nonzero(none_vals) == 1:
            # Must be no iterables that are not strings. Thus, an element-wise
            #  test for iterability and an element-wise test for stringiness
            #  must give matching arrays
            if not all(np.equal(list(map(np.iterable, arglist)), \
                        list(map(lambda e: isinstance(e, str), arglist)))):
                raise ValueError("'None' as parameter invalid with \
                                                        non-str iterables")
            ## end if

            # Parameters okay; replace the None with the appropriate range()
            none_found = True
            none_loc = np.nonzero(none_vals)[0][0]
            arglist[none_loc] = \
                    range(self.num_geoms if none_loc == 0 else self.num_atoms)
        ## end if

        # Return the arguments list and the none-found value
        return arglist

    ## end def _none_subst


    @staticmethod
    def _iter_return(tup, fxn, invalid_error):
        """ Wrapper for exception/|None| output handling of X_iter methods.

        Attempts to pass `tup` as arguments to `fxn`.  If the call is
        successful, returns the value produced. If
        :exc:`~exceptions.IndexError` or
        :exc:`~exceptions.ValueError`
        is raised, indicating an invalid index value, re-raise if
        `invalid_error` is |True| or return |None| if |False|.

        Other exceptions are left uncaught.

        Parameters
        ----------
        tup
            |tuple| --
            Input arguments to unpack as arguments to `fxn`

        fxn
            callable --
            Function/method to call to calculate the return value

        invalid_error
            |bool| --
            Flag for whether to return |None| on an
            :exc:`~exceptions.IndexError` or
            :exc:`~exceptions.ValueError`, or to re-raise.

        Returns
        -------
        val
            |npfloat|_ or |None| --
            Calculated value from fxn(*tup) call, or 'None' value indicating
            IndexError / ValueError

        """

        try:
            val = fxn(*tup)
        except (IndexError, ValueError):
            if invalid_error:
                # Raise the exception if invalid_error indicates
                raise
            else:
                # Otherwise, just return a 'None' value
                return None
            ## end if
        else:
            # Good value; just generate it
            return val
        ## end try

    ## end def _iter_return

## end class OpanXYZ


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable")
