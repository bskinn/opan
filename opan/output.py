#-------------------------------------------------------------------------------
# Name:        output
# Purpose:     Encapsulation of parsing/handling of data from computation
#                textual output.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     16 Feb 2015
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


# Debug constant
_DEBUG = False


class OrcaOutput(object):
    """ Container for parsed textual output generated by ORCA.

    All implemented results that are found in the indicated output are stored
    in the OrcaOutput instance.  If a given quantity was not detectable, it
    is stored as None in the corresponding instance variable.

    The verbose contents of the output file are not generally retained within
    the OrcaOutput instance due to the potential for such to involve a
    tremendously large string. Exceptions include, if present:
        - THERMOCHEMISTRY section

    Instantiation
    -------------
    __init__(output_src, src_type='file')
        Constructor for ORCA output parsing.

    Class Variables
    ---------------
    dict() lookups of re.compile() RegEx patterns:
        p_en    : Energies reported at the end of SCF cycles. Keys:
            EN.SCFFINAL : SCF energy including gCP, D3, etc. corrections
            EN.GCP      : gCP correction
            EN.D3       : D3 correction (D3BJ, at least; unconfirmed with
                            D3ZERO. Likely nonfunctional with DFT-NL)
            EN.SCFOCC   : SCF energy with only COSMO outlying q correction
            EN.OCC      : COSMO outlying q correction
            EN.SCFFINALOCC  : Final SCF energy with COSMO q correction added
        p_thermo    : Quantities extracted from THERMOCHEMISTRY block. Keys:
            THERMO.BLOCK    : Entire THERMOCHEMISTRY block
            THERMO.TEMP     : Simulated temperature (K)
            THERMO.PRESS    : Simulated pressure (atm)
            THERMO.E_EL     : Electronic energy from thermo (Eh; often slightly
                                different than the last EN_SCFFINAL value)
            THERMO.E_ZPE    : Zero-point energy in thermo (Eh)
            THERMO.E_VIB    : Thermal vibrational U correction (Eh)
            THERMO.E_ROT    : Thermal rotational U correction (Eh)
            THERMO.E_TRANS  : Thermal translational U correction (Eh)
            THERMO.H_IG     : Ideal-gas (kB*T) enthalpy contribution (Eh)
            THERMO.TS_EL    : Electronic T*S contribution (Eh)
            THERMO.TS_VIB   : Vibrational T*S contribution (Eh)
            THERMO.TS_TRANS : Translational T*S contribution (Eh)
            THERMO.QROT     : Rotational partition function (unitless)
        p_spincont  : Spin contamination block values. Keys:
            SPINCONT.ACTUAL : Calculated <S**2> expectation value
            SPINCONT.IDEAL  : Ideal <S**2> expectation value for system
            SPINCONT.DEV    : Deviation (calc - ideal)

    str constants for keys of above 'dict struct' variables are also defined,
        as class variables.  For example, to retrieve the Regex pattern for
        locating the vibrational entropy contribution, the following syntax
        can be used:

            OrcaOutput.p_thermo[THERMO.TS_VIB]

    Instance Variables
    ------------------
    src_type : str
        String describing the nature of the source used to create the instance.
    src_src : str
        Descriptor of the location of the source used to create the instance.
    completed : bool
        True if ORCA output reports normal termination, False otherwise.
    converged : bool
        True if SCF converged ANYWHERE IN run. #DOC Update oo.converged
        with any robustifications
    optimized : bool
        True if any OPT converged ANYWHERE in run. Fine for OPT, but ambiguous
        for scans. #DOC Update oo.optimized with any robustifications
    en : dict of lists of np.float64
        Lists of the various energy values from the parsed output. Dict
        keys are those of p_en, above.  Any energy type not found in the
        output is assigned as an empty list.
    thermo : dict of np.float64
        Values from the thermochemistry block of the parsed output. Dict keys
        are those of p_thermo, above.
        #TODO: OrcaOutput.thermo: Test on single-atom case, update above
        #       documentation to reflect outcome
    spincont    : dict of np.float64
        Lists of the various Values from the spin contamination calculations
        in the output, if present. Empty lists if absent. Dict keys are those
        of p_spincont, above.
    thermo_block : str
        Full text of the thermochemistry block, if found.

    Instance Methods
    -------
    en_last : Returns a dict providing the various energy values from the
        last SCF cycle performed in the output. Keys are those of p_en,
        above.  Any energy value not relevant to the parsed output is
        assigned as None.

    Generators
    ----------
    (none)

    """

    # Imports
    import re as _re


    # Various class-level RegEx patterns, collected into dictionaries to
    #  facilitate later iterable data retrieval.
    #
    # IN ALL PATTERNS the group name is the same -- this is to simplify the
    #  parsing process when these patterns are used -- no need to dynamically
    #  fiddle with substituting in custom group names each time! The .replace()
    #  call in each pattern definition saves work if P_GROUP ever needs to be
    #  changed.
    P_GROUP = "val"

    # Patterns for SCF energies, reported at the end of single-point and each
    #  step of geometry optimizations, etc., if not suppressed by %output
    #  settings.

    # String constants for retrieving energy quantities.
    # Prefix is the uppercase of the Regex dictionary name
    class EN(object):
        SCFFINAL = "scffinal"
        GCP = "gcp"
        D3  = "d3"
        SCFOCC = "scfocc"
        OCC = "occ"
        SCFFINALOCC = "scffinalocc"
    ## end class EN

    # Initialize dictionary
    p_en = dict()

    # Final SCF energy, with gCP, D3, corrections included... but NOT COSMO
    #  outlying charge correction.
    p_en.update({ EN.SCFFINAL :
        _re.compile("""
        -\\n                         # Hyphen on preceding line
        FINAL\\ SINGLE\\             # Key text 1
        POINT\\ ENERGY\\             # Key text 2
        [\\ ]+(?P<>[0-9.-]+)         # Energy on same line as key text
        .*\\n-                       # Hyphen starting following line
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })

    # gCP corrections entering into reported FINAL ENERGY values.
    p_en.update({ EN.GCP :
        _re.compile("""
        -\\n                        # Hyphen on preceding line
        gCP\\ correction            # Key text
        [\\ ]+(?P<>[0-9.-]+)        # Energy on same line as key text
        .*\\n-                      # Hyphen starting following line.
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })

    # D3 corrections entering into reported FINAL ENERGY values.
    p_en.update({ EN.D3 :
        _re.compile("""
        -\\n                        # Hyphen on preceding line
        Dispersion\\ correction     # Key text
        [\\ ]+(?P<>[0-9.-]+)        # Energy on same line as key text
        .*\\n-                      # Hyphen starting following line.
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })

    # COSMO SCF energies after the COSMO outlying charge correction BUT BEFORE
    #  any other augmentations to the energy (no D3, gCP, etc.)
    p_en.update({ EN.SCFOCC :
        _re.compile("""
        Total\\ Energy\\ after\\    # Key text 1
        outlying\\ charge\\         # Key text 2
        correction[\\ ]*=           # Key text 3
        [\\ ]+(?P<>[0-9.-]+)        # Energy following key text
        [\\ ]+Eh.*\\n               # 'Eh' units, then newline
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })

    # Patterns for the entire thermochemistry block, as well as the individual
    #  data elements therein.

    # String constants for retrieving energy quantities.
    # Prefix is the uppercase of the Regex dictionary name
    class THERMO(object):
        BLOCK = "block"
        TEMP = "temp"
        PRESS = "press"
        E_EL = "e_el"
        E_ZPE = "e_zpe"
        E_VIB = "e_vib"
        E_ROT = "e_rot"
        E_TRANS = "e_trans"
        H_IG = "h_ig"
        TS_EL = "ts_el"
        TS_VIB = "ts_vib"
        TS_TRANS = "ts_trans"
        QROT = "qrot"
    ## end class THERMO

    # Initialize dictionary
    p_thermo = dict()

    # Whole thermo block just in case; probably not needed for automated
    #  computation, but potentially handy for manual fiddling.
    p_thermo.update({ THERMO.BLOCK :
        _re.compile("""
        (?P<>-+\\n              # Hyphen line
        THERMOCHEMISTRY\\       # Header text
        AT\\ [0-9.]+\\ *K\\n    # Temperature
        -+\\n                   # Hyphen line
        (.|\\n)*)               # Everything until the end
        Timings\\ for\\         # Closing blip 1
        individual\\ modules    # Closing blip 2
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })

    # Individual quantities. Descriptions in pattern definition comments.
    p_thermo.update({ THERMO.TEMP :
        _re.compile("""
        temperature[\\ .]+      # Key text
        (?P<>[0-9.]+)           # Temperature value
        [\\ ]+K                 # in Kelvin
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.PRESS :
        _re.compile("""
        pressure[\\ .]+         # Key text
        (?P<>[0-9.]+)           # Pressure value
        [\\ ]+atm               # in atm
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.E_EL :
        _re.compile("""
        electronic\\ energy     # Key text 1
        [\\ .]+                 # Key text 2
        (?P<>[0-9.-]+)          # Electronic energy value
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.E_ZPE :
        _re.compile("""
        zero\\ point\\ energy   # Key text 1
        [\\ .]+                 # Key text 2
        (?P<>[0-9.-]+)          # ZPE energy value
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.E_VIB :
        _re.compile("""
        thermal\\ vibrational\\ # Key text 1
        correction[\\ .]+       # Key text 2
        (?P<>[0-9.-]+)          # Vibration energy correction value
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.E_ROT :
        _re.compile("""
        thermal\\ rotational\\  # Key text 1
        correction[\\ .]+       # Key text 2
        (?P<>[0-9.-]+)          # Rotation energy correction value
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.E_TRANS :
        _re.compile("""
        thermal\\               # Key text 1
        translational\\         # Key text 2
        correction[\\ .]+       # Key text 3
        (?P<>[0-9.-]+)          # Translation energy correction value
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.H_IG :
        _re.compile("""
        thermal\\ enthalpy\\    # Key text 1
        correction[\\ .]+       # Key text 2
        (?P<>[0-9.-]+)          # Ideal gas enthalpy correction
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.TS_EL :
        _re.compile("""
        electronic\\ entropy\\  # Key text
        [\\ .]+                 # Spacer
        (?P<>[0-9.-]+)          # Electronic entropy contribution
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.TS_VIB :
        _re.compile("""
        vibrational\\ entropy\\ # Key text
        [\\ .]+                 # Spacer
        (?P<>[0-9.-]+)          # Vibrational entropy contribution
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.TS_TRANS :
        _re.compile("""
        translational\\         # Key text 1
        entropy\\               # Key text 2
        [\\ .]+                 # Spacer
        (?P<>[0-9.-]+)          # Translational entropy contribution
        [\\ ]+Eh                # in Hartrees
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_thermo.update({ THERMO.QROT :
        _re.compile("""
        qrot\\ +=\\ +           # Key text
        (?P<>[0-9.]+)           # Rotational partition function
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })


    # Dipole moment pattern
    p_dipmom = _re.compile("""
        -+\\n                             # Hyphen line
        dipole\\ moment.*\\n              # Block label
        -+\\n                             # Hyphen line
        (.*\\n)+?                         # Lazy grab of any lines
        Magnitude\ \(debye\)\ +:\ +       # Line leader
        (?P<>[0-9.]+).*\\n                # Grab the dipole moment
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)


    # Patterns for the spin contamination information
    # String constants for retrieving energy quantities.
    # Prefix is the uppercase of the Regex dictionary name
    class SPINCONT(object):
        ACTUAL = "actual"
        IDEAL = "ideal"
        DEV = "dev"
    ## end class SPINCONTAM

    # Initialize dictionary
    p_spincont = dict()

    # Patterns for the spin contamination information
    p_spincont.update({ SPINCONT.ACTUAL :
        _re.compile("""
        expectation[ ]value[ ]of[ ]<S\*\*2>         # Key text
        [ ]+:[ ]+                                   # Space and separator
        (?P<>[0-9.]+)                               # Grab the value
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_spincont.update({ SPINCONT.IDEAL :
        _re.compile("""
        ideal[ ]value[ ]s\\*\\(s\\+1\\)[ ]          # Key text 1
        for[ ]s=[0-9.]+                             # Key text 2
        [ ]+:[ ]+                                   # Space and separator
        (?P<>[0-9.]+)                               # Grab the value
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })
    p_spincont.update({ SPINCONT.DEV :
        _re.compile("""
        deviation                                   # Key text
        [ ]+:[ ]+                                   # Space and separator
        (?P<>[0-9.]+)                               # Grab the value
        """.replace("P<", "P<" + P_GROUP), _re.I | _re.X)
        })

    #RESUME: Make patterns and constants for the virial block; update docstrings

    ## end class variables

    def __init__(self, output_src, src_type="file"):
        """ Initialize OrcaOutput object.

        Actions performed depend on the indicated 'src_type':

        'file':
            ORCA output is read from a file on disk

        Available data includes:
            SCF energies (incl D3, gCP, COSMO outlying charge corrections)
            Thermochemistry
            Spin expectation values (actual, ideal, and deviation)

        Success indicators include:
            completed : Checks for the 'ORCA TERMINATED NORMALLY' report at the
                    end of the file
            converged : Checks for any occurrence of successful SCF convergence
                    in the file (questionable for anything but single-point
                    calculations)
            optimized : Checks for any occurrence of "OPTIMIZATION HAS
                    CONVERGED" in the file (questionable for anything but
                    a standalone OPT -- i.e., not a mode or internal coordinate
                    scan)

        *** ONLY WORKS ON A VERY SMALL SUBSET OF COMPUTATION TYPES!!! ***

        Parameters
        ----------
        output_src : str
            Depending on the value of 'src_type':
                "file" : Full path to the output file to be parsed.
        src_type : str
            One of the following values indicating the type of output data
            source to be parsed:
                file : ORCA output file on disk

        Raises
        ------
        OUTPUTError : If indicated output is un-parseably malformed in some
                        fashion
        ValueError  : If 'src_type' is invalid
        """
        #TODO: (?) OrcaOutput: Add initialization parameter to indicate which
        # type of run should be expected?

        # Imports
        from .utils import pack_tups
        from .utils import safe_cast as scast
        from .error import OUTPUTError
        import numpy as np

        # Confirm src_type is valid
        if not src_type == "file":
            raise(ValueError("'{0}' is invalid".format(src_type)))
        ##end if

        # Get the output data
        if src_type == "file":
            with open(output_src) as in_f:
                datastr = in_f.read()
            ##end with
        ##end if

        # Check for normal termination (weird values in dicts, etc. would be
        #  diagnostic also, but might as well define this since it's easy).
        self.completed = datastr.find("ORCA TERMINATED NORMALLY") > -1

        # Simple check for single-point SCF convergence
        # TODO: Probably robustify convergence check to opt, and MDCI/MRCI/CAS
        self.converged = datastr.find("SCF CONVERGED AFTER") > -1

        # Simple check for optimization convergence
        # TODO: Probably robustify optimization convergence check for
        #  scans as well as single optimizations.  Multiple job runs promise
        #  to be thoroughly annoying.
        self.optimized = datastr.find("OPTIMIZATION HAS CONVERGED") > -1

        # Store the source information
        self.src_type = src_type
        self.src_src = output_src

        # Initialize the energies dict as empty
        self.en = dict()

        # Populate with the Regex-retrieved values.
        #  If any are not found, this will store as an empty list.
        for (k,p) in self.p_en.iteritems():
            self.en.update({ k :
                    [scast(m.group(self.P_GROUP), np.float_) for m in
                        p.finditer(datastr)] })
        ##next (k,p)

        # Calculate just the outlying charge correction, if COSMO enabled,
        #  and then calculate the SCFFINAL result including the OCC.
        if not self.en[self.EN.SCFOCC] == []:
            self.en.update({ self.EN.OCC :
                    [t[0] - (t[1] - t[2] - t[3]) for t in
                    pack_tups(
                        self.en[self.EN.SCFOCC],
                        self.en[self.EN.SCFFINAL],
                        self.en[self.EN.D3] if self.en[self.EN.D3] != []
                                                                else 0,
                        self.en[self.EN.GCP] if self.en[self.EN.GCP] != []
                                                                else 0
                              )
                    ]       })
            self.en.update({ self.EN.SCFFINALOCC :
                    [t[0] + t[1] for t in
                    pack_tups( # Could use zip() here, probably
                        self.en[self.EN.SCFFINAL],
                        self.en[self.EN.OCC]
                              )
                    ]       })
        ##end if

        # Now collect the thermo quantities
        # Just store the whole thermo block
        try:
            self.thermo_block = \
                    self.p_thermo[self.THERMO.BLOCK].search(datastr).group()
        except AttributeError:
            # Block not found; store as None
            self.thermo_block = None
        else:
            # Only store the block details if the block is actually found!

            # Initialize the empty dictionary for the numericals
            self.thermo = dict()

            # Iterate to pull the individual values
            for (k,p) in self.p_thermo.iteritems():
                if k != self.THERMO.BLOCK:
                    try:
                        self.thermo.update({ k : \
                                    scast(p.search(datastr) \
                                            .group(self.P_GROUP), np.float_) })
                    except AttributeError:
                        # Value not found, probably due to monoatomic freq calc
                        #  to autogenerate, e.g., enthalpy calculation
                        # Add as just None'
                        self.thermo.update({ k: None })
                    ## end try
                ## end if
            ## next (k,p)
        ## end try

        #TODO: (?) OrcaOutput: Pull the final geometry and atom masses. Would
        #  be nice not to require a Hessian calculation in order to have this
        #  info available.
        #  Masses and/or geometries may not always be in the output file,
        #  depending on the %output settings.  Also have to address possible
        #  multiples of the coordinates in, e.g., scans.

        # Pull all dipole moments
        self.dipmoms = []
        for m in OrcaOutput.p_dipmom.finditer(datastr):
            self.dipmoms.append(scast(m.group(OrcaOutput.P_GROUP), np.float_))
        ## next m

        # Initialize the spin contamination dict as empty
        self.spincont = dict()

        # Populate with the Regex-retrieved values.
        #  If any are not found, this will store as an empty list.
        for (k,p) in self.p_spincont.iteritems():
            self.spincont.update({ k :
                    [scast(m.group(self.P_GROUP), np.float_) for m in
                        p.finditer(datastr)] })
        ##next (k,p)

        #RESUME: Pull the virial block info (may be absent)

    ## end def __init__


    def en_last(self):
        """ Report the energies from the last SCF present in the output.

        Returns a dict providing the various energy values from the
        last SCF cycle performed in the output. Keys are those of p_en; see
        OrcaOutput.__doc__.  Any energy value not relevant to the parsed
        output is assigned as None.

        Parameters
        ----------
        (none)

        Returns
        -------
        last_ens : dictionary of np.float64
            Energies from the last SCF present in the output. Reports all
            energy values implemented. Keys are identical to, e.g., the RegEx
            pattern dicts used for searching the output; see the class
            docstring.

        Raises
        ------
        (none)
        """

        # Initialize the return dict
        last_ens = dict()

        # Iterate and store
        for (k,l) in self.en.items():
            last_ens.update({ k : l[-1] if l != [] else None })
        ##next (k,l)

        # Should be ready to return?
        return last_ens

    ## end def en_last


if __name__ == '__main__':
    print("Module not executable")
