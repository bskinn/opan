#-------------------------------------------------------------------------------
# Name:        error
# Purpose:     Definitions of all custom Errors for the OpenAnharmonic
#                package.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     1 Oct 2014
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

# Module-level imports


class OPANError(Exception):
    """Base class for custom Errors defined for OpenAnharmonic

    OPANError is an abstract superclass of any custom errors built under the
    OpenAnharmonic umbrella. It defines all common methods shared among the
    various subtype error classes, such that the only contents that must be
    declared by a subclass (as of this comment) is the 'typecodes' class-level
    dictionary, which maps string typecode keys to *unique* typecode int values.

    Instantiation
    -------------
    tc  : str
        String representation of typecode to be associated with the OPANError
        subclass instance. *Must* be a valid key of the 'typecodes' class
        variable defined for the relevant subtype.
    msg : str
        Explanation of the nature of the error being reported.
    src : str
        Further detail of the code source of the error behavior.

    Class Variables
    ---------------
    typecodes: dict
        (Defined separately by each subclass of OPANError.) Dictionary of the
        allowed typecodes for the relevant OPANError subclass.

    Instance Variables
    ------------------
    msg : str
        Explanation of the nature of the error being reported.
    src : str
        Further detail of the code source of the error behavior.
    subclass_name   : str
        String representation of the name of the OPANError subclass of the
        current OPANError instance.
    tc  : int
        Numeric value of the typecode associated with the instance, as defined
        in the 'typecodes' dict of the subclass

    Class Methods
    -------
    _typecode_str -- Returns string representation of a numerical typecode
        value. Meant primarily for introspection.

    """


    def __init__(self, tc, msg, src):
        """ Uniform constructor for subclasses of abstract OPANError.

        Parameters
        ----------
        tc   :  typecode for error
                    (see subclass docstrings for allowed typecodes)
        msg  :  explanation of the error
        src  :  source of the problematic data

        Returns
        -------
        (none)

        Raises
        ------
        KeyError : Invalid typecode provided in 'tc'
        AttributeError : If ORCAError subclass has not defined the class
                            variable 'typecodes'
        NotImplementedError : Upon attempt to instantiate abstract ORCAError
            base class.

        """

        # Import(s)
        import re

        # Check for and complain at instantiation of base class. Will need
        #  to update this every time a new subclass of error is created.
        if not (
                    isinstance(self, XYZError) or \
                    isinstance(self, GRADError) or \
                    isinstance(self, HESSError) or \
                    isinstance(self, OUTPUTError) or \
                    isinstance(self, SYMMError) or \
                    isinstance(self, REPOError) or \
                    isinstance(self, ANHARMError)
                ):
            raise(NotImplementedError("OPANError base class is abstract."))
        ## end if

        # Quick RegEx to extract the name of the subclass.
        self.subclass_name = re.search(self.__module__ + "\\.(?P<cls>\w+)'", \
                    str(self.__class__), re.I).group("cls")

        # Check for valid typecode and throw a more descriptive error if
        #  invalid.
        if not tc in self.typecodes:
            raise(KeyError("Invalid " + self.subclass_name + \
                                                " typecode: " + str(tc)))
        ## end if

        # Store error content
        self.tc = tc
        self.msg = msg
        self.src = src


    def __str__(self):
        """ String representation of an instance of an ORCAError subclass.

        Implemented primarily so that the error stack handling of the Python
        interpreter will provide useful information to the user.

        Return value is constructed as follows:
            (typecode string) Error message: Error source

        Parameters
        ----------
        (none)

        Returns
        -------
        retstr : str
            String representation of the instance.

        Raises
        ------
        (none)

        """

        # Store and return the descriptive string
        retstr = "(" + str(self.tc) + ") " + str(self.msg) + \
                (": " + str(self.src) if str(self.src) != "" else "")
        return retstr

    ## end def __str__


    class __metaclass__(type):
        # Enable iteration over the typecodes
        def __iter__(self):
            for item in self.__dict__:
                if item == self.__dict__[item]:
                    yield item
                ## end if
            ## next item
        ##end def __iter__
    ## end class __metaclass__


    @classmethod
    def _typecode_str(self, tc):
        """Return string representation of provided typecode number.

        Helper function for introspection lookup of the typecode string of the
        typecode number input as 'tc'.  Raises ValueError if the typecode value
        is not found in the .typecodes dict of the appropriate subclass.

        POTENTIALLY OBSOLETE!

        Parameters
        ----------
        tc  : int
            Typecode enumeration constant to be converted to its string
            representation

        Returns
        -------
        tc_str  : string
            Representation of the indicated typecode for the class of which
            'self' is a member.

        Raises
        ------
        ValueError  : If an invalid typecode value is passed.

        """

        # List comprehension(?) lookup of the indicated typecode
        try:
            tc_str = str([k for (k,v) in self.typecodes.items() if v == tc][0])
        except IndexError:
            raise(ValueError("Invalid typecode value."))
        ##end try

        return tc_str

## end class ORCAError


class XYZError(OPANError):
    """Error relating to parsing of or calculation from XYZ data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            xyzfile     :  inconsistent geometry in an OpenBabel XYZ
                            file (XYZ or TRJ from ORCA)
            overwrite   :  object already initialized (overwrite not supported)
            dihed       :  dihedral angle calculation requested for a set
                            of atoms containing an insufficiently
                            nonlinear trio of atoms
            nonprl      :  insufficient non-parallel character in some
                            manner of calculation

    """

    # Typecodes as class-level variables, collected into frozenset
    xyzfile = 'xyzfile'
    overwrite = 'overwrite'
    dihed = 'dihed'
    nonprl = 'nonprl'

    typecodes = frozenset([
            xyzfile,
            overwrite,
            dihed,
            nonprl
            ])

## end class XYZError


class GRADError(OPANError):
    """Error relating to parsing of or calculation from gradient data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            numats     : Invalid number of atoms specification,
                            or specification not found
            en         : Energy specification not found
            gradblock  : Invalid gradient block or gradient block
                            not found
            geomblock  : Invalid geometry block or geometry block not found

    """

    # Typecodes as class-level variables, collected into frozenset
    numats = 'numats'
    en = 'en'
    gradblock = 'gradblock'
    geomblock = 'geomblock'

    typecodes = frozenset([
            numats,
            en,
            gradblock,
            geomblock
            ])

## end class ENGRADError


class OUTPUTError(OPANError):
    """Error relating to parsing of or calculation from output data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            {add here}

    """

    # Typecodes as class-level variables, collected into frozenset
    typecodes = frozenset([

            ])

## end class OUTPUTError


class HESSError(OPANError):
    """Error relating to parsing of or calculation from Hessian data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            at_block    : Malformed or missing atom specification block
            energy      : Malformed or missing energy value
            freq_block  : Malformed or missing frequencies block
            hess_block  : Malformed or missing Hessian block
            modes_block : Malformed or missing normal modes block
            temp        : Malformed or missing 'actual temperature' value

            dipder_block: Malformed dipole derivatives block
            eigval_block: Malformed mass-wt Hessian eigenvalues block
            eigvec_block: Malformed mass-wt Hessian eigenvectors block
            ir_block    : Malformed IR spectrum block
            job_block   : Malformed job list block
            polder_block: Malformed polarizability derivatives block
            raman_block : Malformed Raman spectrum block


    """

    # Typecodes as class-level variables, collected into frozenset
    at_block = "at_block"
    hess_block = "hess_block"
    freq_block = "freq_block"
    modes_block = "modes_block"
    dipder_block = "dipder_block"
    ir_block = "ir_block"
    polder_block = "polder_block"
    raman_block = "raman_block"
    job_block = "job_block"
    eigval_block = "eigval_block"
    eigvec_block = "eigvec_block"
    energy = "energy"
    temp = "temp"
    typecodes = frozenset([
            at_block,
            hess_block,
            freq_block,
            modes_block,
            dipder_block,
            ir_block,
            polder_block,
            raman_block,
            job_block,
            eigval_block,
            eigvec_block,
            energy,
            temp
            ])

## end class HESSError


class SYMMError(OPANError):
    """Error relating to determination of molecular symmetry.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            notfound    : Symmetry element expected but not found.

    """

    # Typecodes as class-level variables, collected into frozenset
    notfound = 'notfound'

    typecodes = frozenset([
            notfound
            ])

## end class SYMMError


class REPOError(OPANError):
    """Error relating to HDF5 repository interactions.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            status      : HDF5 repo in improper status for requested operation.
            group       : Problem with a group in linked HDF5 file
            data        : Problem with a dataset in linked HDF5 file

    """

    # Typecodes as class-level variables, collected into frozenset
    status = 'status'
    data = 'data'
    group = 'group'

    typecodes = frozenset([
            data,
            status,
            group
            ])

## end class REPOError


class ANHARMError(OPANError):
    """Error relating to OPAN_ANHARM class actions.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes: dict
            repo        : OPAN_REPO conflict -- no repo bound when assignment
                            attempted, or attempt made to bind new repo when
                            one already bound
            status      : OPAN_ANHARM internal variables in inappropriate
                            status for the requested operation

    """

    # Typecodes as class-level variables, collected into frozenset
    repo = 'repo'
    status = 'status'

    typecodes = frozenset([
            repo,
            status
            ])

## end class ANHARMError


if __name__ == '__main__':
    print("Module not executable")
