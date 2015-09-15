#-------------------------------------------------------------------------------
# Name:        hess
# Purpose:     Encapsulation of parsing/handling of data from Hessian files.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     14 Nov 2014
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

# Imports


# Debug constant
_DEBUG = False


class ORCA_HESS(object):
    """ Container for HESS data generated by ORCA.

    Information contained includes the Hessian matrix, the number of atoms,
    the atomic symbols, the atomic weights, and the geometry, as reported in
    the .hess file.  The precision of the geometry is less than that reported
    in an .xyz file, and thus should NOT be used for generation of subsequent
    computations.

    'N' in the below documentation refers to the number of atoms present in the
    geometry contained within the ENGRAD.

    Constructor may need to be adapted at some point in the future to handle
    construction from a custom storage container. For now, the file-based
    retrieval is the only implemented mode. If h5py turns out to be well
    suited for use as a repository format, such modification may not be
    necessary.

    Units of the Hessian are Hartrees per Bohr^2 (Eh/B^2)
    Frequencies are in cyc/cm (standard wavenumbers)


    Instantiation
    -------------
    __init__(HESS_path)
        Constructor for an ORCA_HESS, drawing data from a HESS file on disk.

    Class Variables
    ---------------
    p_at_block      : re.compile() pattern
        Regex for the entire block of atom ID & weight, and geom data.
    p_at_line       : re.compile() pattern
        Regex for individual lines within the atom specification block
    p_energy        : re.compile() pattern
        Regex for the energy reported in the .hess file
    p_freq_block    : re.compile() pattern
        Regex for the entire vibrational frequencies block (cm**-1 units)
    p_freq_line     : re.compile() pattern
        Regex for each line of the frequencies block
    p_hess_block    : re.compile() pattern
        Regex for entire Hessian block
    p_hess_sec      : re.compile() pattern
        Regex for a full-height, 3- or 6-column section of the Hessian
    p_hess_line     : re.compile() pattern
        Regex for a single line within a Hessian section
    p_ir_block      : re.compile() pattern
        Regex for the full IR spectrum block
    p_ir_line       : re.compile() pattern
        Regex for individual IR spectrum lines
    p_modes_block   : re.compile() pattern
        Regex for entire modes block
    p_modes_sec     : re.compile() pattern
        Regex for a full-height, 3- or 6-column section of the modes block
    p_modes_line    : re.compile() pattern
        Regex for a single line within a modes block section
    p_dipder_block  : re.compile() pattern
        Regex for the dipole derivatives block
    p_dipder_line   : re.compile() pattern
        Regex for individual lines in the dipole derivatives block
    p_temp          : re.compile() pattern
        Regex for the 'actual temperature' field (may be meaningless?)

    Instance Variables
    ------------------
    atom_masses     : N x 1 np.float_
        Column vector of atom masses as reported in HESS
    atom_syms       : N x 1 np.str
        Column vector of uppercase atomic symbols
    dipders         : 3N x 3 np.float_
        Matrix of dipole derivatives
    energy          : float
        Energy reported in the Hessian file
    freqs           : 3N x 1 np.float_
        Vibrational frequencies (cm**-1) as reported in the Hessian file
    geom            : 3N x 1 np.float_
        Column vector of geometry [x1, y1, z1, x2, y2, ...]
    hess            : 3N x 3N np.float_
        Cartesian Hessian matrix
    HESS_path       : str
        Complete path/filename from which the Hessian data was retrieved
    initialized     : bool
        Flag for whether self has been initialized--possibly obsolete
    in_str          : str
        Complete contents of the imported HESS file
    modes           : 3N x 3N np.float_
        Rotation- and translation-purified vibrational normal modes,
        with each mode (column vector) individually normalized by ORCA.
    num_ats         : int
        Number of atoms in the system
    temp            : float
        "Actual temperature" reported in the .hess file. May be meaningless.

    Methods
    -------
    check_geom(coords, atoms[, tol])
        Checks vectors of atom coordinates and identities for consistency
        with the geometry and atom identities stored within the instance of
        ORCA_HESS.

    Generators
    ----------
    (none)
    """

    # Imports
    import re
    from .const import DEF


    # Various class-level Regex patterns.  Currently only retrieves the
    #  atom list and the Hessian itself. No other information in the .hess file
    #  is actually explicitly needed at present (will be recomputed internally)
    #  and so there's little point to spending time coding imports for it.

    # Atoms list, including atomic weights. Assumes no scientific notation
    #  will be used in the coordinates.
    p_at_block = re.compile("""
    \\#.*\\n                        # Line prior to $atoms is a blank comment
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
    """, re.I | re.M | re.X)

    # Pulling the individual atoms, weights, and coordinates from each atom line
    p_at_line = re.compile("""
    ^[ \\t]+                        # Whitespace to start the line
    (?P<el>([a-z]+|[0-9]+))         # Atomic number or element symbol
    [ \\t]+(?P<mass>[0-9.]+)        # Whitespace and atomic mass
    [ \\t]+(?P<c1>[0-9.-]+)         # Whitespace and first coordinate
    [ \\t]+(?P<c2>[0-9.-]+)         # Whitespace and second coordinate
    [ \\t]+(?P<c3>[0-9.-]+)         # Whitespace and third coordinate
    """, re.I | re.M | re.X)

    # Entire Hessian data block
    p_hess_block = re.compile("""
    \\$hessian.*\\n                 # Marker for Hessian block
    (?P<dim>[0-9]+).*\\n            # Dimensionality of Hessian (3N x 3N)
    (?P<block>                      # Group for the subsequent block of lines
        (                           # Group for single line definition
            ([ \\t]+[0-9.-]+)+      # Some number of whitespace-separated nums
            .*\\n                   # Plus whatever to end of line
        )+                          # Whatever number of single lines
    )                               # Enclose the whole batch of lines
    """, re.I | re.X)

    # Sections of the Hessian data block
    p_hess_sec = re.compile("""
    ([ \\t]+[0-9]+)+[ \\t]*\\n      # Column header line
    (                               # Open the group for the sub-block lines
        [ \\t]+[0-9]+               # Row header
        (                           # Open the group defining a single element
            [ \\t]+[-]?             # Whitespace and optional hyphen
            [0-9]+\\.[0-9]+         # One or more digits, decimal, more digits
        )+                          # Some number of sub-columns
        [ \\t]*\\n                  # Whitespace to EOL
    )+                              # Some number of suitable lines
    """, re.I | re.X)

    # Pulling Hessian lines from the sections, with elements in groups
    p_hess_line = re.compile("""
    ^[ \\t]*                        # Optional whitespace to start each line
    (?P<row>[0-9]+)                         # Row header
    [ \\t]+(?P<e0>[0-9-]+\\.[0-9]+)         # 1st element
    [ \\t]+(?P<e1>[0-9-]+\\.[0-9]+)         # 2nd element
    [ \\t]+(?P<e2>[0-9-]+\\.[0-9]+)         # 3rd element
    ([ \\t]+(?P<e3>[0-9-]+\\.[0-9-]+))?     # 4th element (possibly absent)
    ([ \\t]+(?P<e4>[0-9-]+\\.[0-9-]+))?     # 5th element (possibly absent)
    ([ \\t]+(?P<e5>[0-9-]+\\.[0-9-]+))?     # 6th element (possibly absent)
    .*$                             # Whatever to end of line
    """, re.I | re.M | re.X)

    # Reported energy
    p_energy = re.compile("""
    \\$act_energy[ ]*\\n                    # Label for the block
    [ ]+(?P<en>[0-9.-]+)[ ]*\\n             # Energy value
    """, re.I | re.X)

    # Frequencies block
    p_freq_block = re.compile("""
    \\$vibrational_frequencies[ ]*\\n       # Block label
    [ ]*(?P<num>[0-9]+)[ ]*\\n              #--> Number of frequencies
    (?P<block>                              # Open capture group for block
        ([ ]+[0-9]+[ ]+[0-9.-]+[ ]*\\n)+    #--> Block contents
    )                                       # Close capture group
    """, re.I | re.X)

    p_freq_line = re.compile("""
    ^[ ]+(?P<id>[0-9]+)                     #--> ID for the frequency
    [ ]+(?P<freq>[0-9.-]+)                  #--> Freq value in cm**-1
    [ ]*$                                   # To EOL
    """, re.I | re.M | re.X)


    # Entire modes data block
    p_modes_block = re.compile("""
    \\$normal_modes.*\\n            # Marker for modes block
    (?P<dim>[0-9]+)[ ]+             # Dimensionality of modes block (3N x 3N)
    (?P<dim2>[0-9]+).*\\n           #  (Second dimension value)
    (?P<block>                      # Group for the subsequent block of lines
        (                           # Group for single line definition
            ([ ]+[0-9.-]+)+         # Some number of whitespace-separated nums
            .*\\n                   # Plus whatever to end of line
        )+                          # Whatever number of single lines
    )                               # Enclose the whole batch of lines
    """, re.I | re.X)

    # Sections of the modes data block
    p_modes_sec = re.compile("""
    ([ ]+[0-9]+)+[ ]*\\n            # Column header line
    (                               # Open the group for the sub-block lines
        [ ]+[0-9]+                  # Row header
        (                           # Open the group defining a single element
            [ ]+[-]?                # Whitespace and optional hyphen
            [0-9]+\\.[0-9]+         # One or more digits, decimal, more digits
        )+                          # Some number of sub-columns
        [ ]*\\n                     # Whitespace to EOL
    )+                              # Some number of suitable lines
    """, re.I | re.X)

    # Pulling modes lines from the sections, with elements in groups
    #  THE USE of the '[0-9-]+\\.[0-9]+' construction here, and in its variants
    #  above/below, guarantees a floating-point value is found. Otherwise, the
    #  Regex retrieves on into subsequent sections because the header rows
    #  parse just fine for a '[ ]+[0-9.-]' pattern.
    p_modes_line = re.compile("""
    ^[ ]*                           # Optional whitespace to start each line
    (?P<row>[0-9]+)                         # Row header
    [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            # 1st element
    [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            # 2nd element
    [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            # 3rd element
    ([ ]+(?P<e3>[0-9-]+\\.[0-9-]+))?        # 4th element (possibly absent)
    ([ ]+(?P<e4>[0-9-]+\\.[0-9-]+))?        # 5th element (possibly absent)
    ([ ]+(?P<e5>[0-9-]+\\.[0-9-]+))?        # 6th element (possibly absent)
    .*$                             # Whatever to end of line
    """, re.I | re.M | re.X)

    # "Actual temperature"
    p_temp = re.compile("""
    \\$actual_temperature[ ]*\\n            # Marker for value
    [ ]*(?P<temp>[0-9.]+)[ ]*\\n            # Value
    """, re.I | re.X)

    # Dipole derivatives block
    p_dipder_block = re.compile("""
    \\$dipole_derivatives[ ]*\\n            # Marker for block
    (?P<dim>[0-9]+)[ ]*\\n                  #--> Dimension of block (rows)
    (?P<block>                              #--> Catches entire block
        (([ ]+[0-9.-]+)+[ ]*\\n)+           # Rows of data
    )                                       # End block capture
    """, re.I | re.X)

    # Dipole derivatives individual line
    p_dipder_line = re.compile("""
    ^                                       # Line start
    [ ]+(?P<e0>[0-9-]+\\.[0-9]+)            # 1st element
    [ ]+(?P<e1>[0-9-]+\\.[0-9]+)            # 2nd element
    [ ]+(?P<e2>[0-9-]+\\.[0-9]+)            # 3rd element
    [ ]*$
    """, re.I | re.M | re.X)


    def __init__(self, HESS_path):
        """ Initialize ORCA_HESS Hessian object from .hess file

        Searches indicated file for geometry and Hessian block. Currently
        does *not* retrieve any of the other information stored in a .hess
        file.

        In future, will likely require extension to handle import from some
        manner of custom object, in order to implement saving/loading of
        anharmonic computations.


        Parameters
        ----------
        HESS_path : string
            Complete path to the .hess file to be read.

        Raises
        ------
        HESSError   : If indicated Hessian file is malformed in some fashion
        KeyError    : If invalid atomic symbol appears in .hess file
        IOError     : If the indicated file does not exist or cannot be read

        Local Methods
        -------------
        parse_multiblock #DOC:Complete parse_multiblock docstring entry
        """

        # Local method(s)
        def parse_multiblock(hesstext, p_block, p_sec, p_line, num_ats, \
                                                            blockname, tc):
            """ #DOC: parse_multiblock docstring
            """

            # Imports
            import numpy as np
            from .utils import safe_cast as scast

            # Pull the block
            m_block = p_block.search(hesstext)

            # Confirm the anticipated matrix size matches 3*num_ats
            if not int(m_block.group("dim")) == 3*num_ats:
                raise(HESSError(tc, blockname + " dimension " + \
                        "specification mismatched with geometry", \
                        "HESS File: " + HESS_path))
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
                        val = scast(m_line.group('e' + str(i)), np.float_)

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
                raise(HESSError(tc, "Insufficient number of " + blockname + \
                        " sections found", \
                        "HESS File: " + HESS_path))
            ## end if

            # Additional cross-check on number of rows imported
            if (rows_counter % (3*num_ats)) != 0:
                # Not the right number of rows; complain
                raise(HESSError(tc, \
                            blockname + " row count mismatch", \
                            "HESS File: " + HESS_path))
            ## end if

            # Return the working matrix as a matrix
            workmtx = np.matrix(workmtx)
            return workmtx

        ## end def parse_multiblock

        # Imports
        import numpy as np
        from .error import HESSError
        from .utils import safe_cast as scast
        from .const import atomNum, atomSym, PRM

        # Set the initialization flag (possibly unnecessary...)
        self.initialized = False

        # Open file, read contents, close stream
        # No particular exception handling; that will be the responsibility
        #  of the calling routine.
        # May want to add a check for text file format; but, really, a check
        #  for expected-information-not-found will likely cover such a case.
        with open(HESS_path,'rU') as in_fl:
            self.in_str = in_fl.read()
        ## end with

        # Store the path used to retrieve the data
        self.HESS_path = HESS_path

        # Check to ensure all required data blocks are found
        if not ORCA_HESS.p_at_block.search(self.in_str):
            raise(HESSError(HESSError.at_block,
                    "Atom specification block not found",
                    "HESS File: " + HESS_path))
        ## end if
        if not ORCA_HESS.p_hess_block.search(self.in_str):
            raise(HESSError(HESSError.hess_block,
                    "Hessian block not found",
                    "HESS File: " + HESS_path))
        ## end if
        if not ORCA_HESS.p_freq_block.search(self.in_str):
            raise(HESSError(HESSError.freq_block,
                    "Frequencies block (cm**-1 units) not found",
                    "HESS File: " + HESS_path))
        ## end if
        if not ORCA_HESS.p_modes_block.search(self.in_str):
            raise(HESSError(HESSError.modes_block,
                    "Normal modes block not found",
                    "HESS File: " + HESS_path))
        ## end if


        #=== Geometry spec block ===#
        # Bring in the number of atoms
        self.num_ats = np.int_( \
                ORCA_HESS.p_at_block.search(self.in_str).group("num"))

        # Initialize the vectors of atomic symbols, masses, and the
        #  geometry
        self.atom_syms = []
        self.atom_masses = []
        self.geom = []

        # Parse the list of atoms
        for m in ORCA_HESS.p_at_line.finditer( \
                ORCA_HESS.p_at_block.search(self.in_str).group("block")):
            # Parse the element symbol or atomic number
            try:
                # See if it casts as an int
                num = scast(m.group('el'), np.int_)
            except ValueError, TypeError:
                # Nope, must be letters. Check if valid by poking it into the
                #  dict. If not valid, should raise another error.
                num = atomNum[str.upper(m.group('el'))]
            ## end try

            # If this point is reached, either it's cast ok to an int,
            #  or the string has been passed through the atomNum dict and
            #  is a valid element symbol. SO, convert back to atomSym and
            #  store - this will also check for an invalid int entry
            self.atom_syms.append(atomSym[num])

            # Parse the atomic weight; should be a simple append?
            self.atom_masses.append(scast(m.group("mass"), np.float_))

            # Build geometry as a 1-D list
            for gp in ["c1", "c2", "c3"]:
                self.geom.append(scast(m.group(gp), np.float_))
            ## next gp

        # Double-check that the number of atoms retrieved matches the
        #  number indicated in the HESS file; geometry size also.
        if not len(self.atom_syms) == self.num_ats:
            raise(HESSError(HESSError.at_block, "Atomic symbol dimension " + \
                    "mismatch with HESS atom specification", \
                    "HESS File: " + HESS_path))
        ## end if
        if not len(self.atom_masses) == self.num_ats:
            raise(HESSError(HESSError.at_block, "Atomic mass dimension " + \
                    "mismatch with HESS atom specification", \
                    "HESS File: " + HESS_path))
        ## end if
        if not len(self.geom) == 3 * self.num_ats:
            raise(HESSError(HESSError.at_block, "Geometry dimension " + \
                    "mismatch with HESS atom specification", \
                    "HESS File: " + HESS_path))
        ## end if

        # Convert instance variables to numpy storage forms
        self.atom_syms = np.matrix(np.vstack(self.atom_syms))
        self.atom_masses = np.matrix(np.vstack(np.array(self.atom_masses, \
                                dtype=np.float_)))
        self.geom = np.matrix(np.vstack(np.array(self.geom, dtype=np.float_)))

        #=== Hessian and modes ===#
        # Pull the Hessian
        self.hess = parse_multiblock(self.in_str, self.p_hess_block, \
                self.p_hess_sec, self.p_hess_line, self.num_ats, \
                "Hessian", HESSError.hess_block)

        # Pull the modes
        self.modes = parse_multiblock(self.in_str, self.p_modes_block, \
                self.p_modes_sec, self.p_modes_line, self.num_ats, \
                "modes", HESSError.modes_block)

        # Extra check of 'dim' vs 'dim2' on modes
        if scast(self.p_modes_block.search(self.in_str) \
                                    .group("dim"), np.int_) != \
                scast(self.p_modes_block.search(self.in_str) \
                                    .group("dim2"), np.int_):
            raise(HESSError(HESSError.modes_block, \
                    "Normal mode block dimension specification mismatch",
                    "HESS File: " + self.HESS_path))
        ## end if

        #=== Frequencies ===#
        # Check that number of frequencies indicated in the block matches
        #  that expected from the number of atoms
        if 3*self.num_ats != \
                np.int_(self.p_freq_block.search(self.in_str).group("num")):
            raise(HESSError(HESSError.freq_block, \
                    "Count in frequencies block != 3 * number of atoms", \
                    "HESS File: " + self.HESS_path))
        ## end if

        # Retrieve the frequencies (this generates a row vector)
        self.freqs = np.matrix( \
                [np.float_(m.group("freq")) for m in \
                self.p_freq_line.finditer( \
                self.p_freq_block.search(self.in_str).group("block"))])

        # Proofread for proper size (still a row vector)
        if not self.freqs.shape[1] == 3*self.num_ats:
            raise(HESSError(HESSError.freq_block, \
                    "Number of frequencies != 3 * number of atoms", \
                    "HESS File: " + self.HESS_path))
        ## end if

        # Transpose for column vector storage
        self.freqs = self.freqs.transpose()


        #=== Pull the single values that should always be present ===#
        # Store the reported energy
        self.energy = scast(self.p_energy.search(self.in_str).group("en"), \
                                                                    np.float_)

        # Store the reported 'actual temperature'
        self.temp = scast(self.p_temp.search(self.in_str).group("temp"), \
                                                                    np.float_)

        #=== Dipole derivatives ===#
        # Check if block found. Store None if not; otherwise import
        if not ORCA_HESS.p_dipder_block.search(self.in_str):
            self.dipders = None
        else:
            # Check that number of derivatives rows indicated in the block
            #  matches that expected from the number of atoms
            if 3*self.num_ats != np.int_(self.p_dipder_block \
                                        .search(self.in_str).group("dim")):
                raise(HESSError(HESSError.dipder_block, \
                        "Count in dipole derivatives block != 3 * # of atoms", \
                        "HESS File: " + self.HESS_path))
            ## end if

            # Retrieve the derivatives
            self.dipders = np.matrix( \
                    [[np.float_(m.group("e" + str(i))) for i in range(3)] \
                        for m in self.p_dipder_line.finditer( \
                            self.p_dipder_block.search(self.in_str) \
                                                            .group("block") \
                                                            )])

            # Proofread for proper size. Don't have to proofread the width of
            #  three, since any row not containing three numerical values will
            #  result in the block getting truncated.
            if not self.dipders.shape[0] == 3*self.num_ats:
                raise(HESSError(HESSError.dipder_block, \
                        "Number of dipole derivative rows != " + \
                                                    "3 * number of atoms", \
                        "HESS File: " + self.HESS_path))
            ## end if

            # If max-absolute element is too big, overwrite with None
            if np.max(np.abs(self.dipders)) > PRM.Max_Sane_DipDer:
                self.dipders = None
            ## end if
        ## end if


        # Set initialization flag; probably unnecessary?
        self.initialized = True

    ## end def __init__


    def check_geom(self, coords, atoms, tol=DEF.HESS_Coord_Match_Tol):
        """ Check for consistency of HESS geometry with input coords/atoms.

        HESS cartesian coordinates are considered consistent with the input
            coords if each component matches to within 'tol' (default value
            specified by orca_const.DEF.HESS_Coord_Match_Tol).  If coords or
            atoms vectors are passed that are of different length than those
            stored in the ORCA_ENGRAD instance, a False value is returned.

        The coords vector must be three times the length of the atoms vector
            or a ValueError is raised.

        Parameters
        ----------
        coords : 3N x 1 float
            Vector of stacked 'lab-frame' Cartesian coordinates
        atoms  : N x 1 string or int
            Vector of atom symbols or atomic numbers
        tol    : float, optional
            Tolerance for acceptable deviation of each geometry coordinate
            from that in the ORCA_HESS instance to still be considered
            matching

        Returns
        -------
        match  : bool
            Whether input coords and atoms match those in the ORCA_ENGRAD
            instance (True) or not (False)
        fail_type  : string
            If match == False, a string description code for the reason
            for the failed match:
                coord_dim_mismatch  : Mismatch in coordinate vector sizes
                atom_dim_mismatch   : Mismatch in atom symbol vector sizes
                coord_mismatch      : Mismatch in one or more coordinates
                atom_mismatch       : Mismatch in one or more atoms
                #DOC: Propagate info when mismatch code converted to Enum
        fail_loc   : 3N x 1 bool or N x 1 bool
            np.matrix() column vector indicating positions of mismatch in
            either coords or atoms, depending on the value of fail_type.
            True elements indicate corresponding *MATCHING* values; False
            elements mark *MISMATCHES*.

        Raises
        ------
        ValueError : If len(coords) != 3 * len(atoms)
        """

        # Import(s)
        from .utils import check_geom as ucg

        # Wrapper call
        result = ucg(self.geom, self.atom_syms, coords, atoms, tol=tol)

        # Return result
        return result

    ## end def check_geom


if __name__ == '__main__':
    print("Module not executable.")


