#-------------------------------------------------------------------------------
# Name:        error
# Purpose:     Definitions of all custom errors for the OpenAnharmonic
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



"""Custom errors for OpenAnharmonic

Error classes are subclassed from :exc:`Exception` via an abstract superclass,
:class:`~opan.error.OPANError`, which defines several common features:

 * Storage of a 'typecode' and a 'source' string along with the error message
   to allow passing of more fine-grained information to the exception stack
 * Implementation of a :attr:`~opan.error.OPANError.__metaclass__` on
   :class:`OPANError` enabling typecode validity checking with `is`
 * Re-implementation of :func:`~OPANError.__str__` to enhance the usefulness
   of stack messages when one of these errors is raised

:class:`OPANError` **Subclasses**

    :class:`ANHARMError` -- Raised as a result of
    :class:`~opan.anharm.OPAN_ANHARM` actions

    :class:`GRADError` -- Raised during parsing of or calculations using
    gradient data

    :class:`HESSError` -- Raised during parsing of or calculations using
    Hessian data

    :class:`INERTIAError` -- Raised by :mod:`opan.utils.inertia`
    submodule functions

    :class:`OUTPUTError` -- Raised during parsing of or calculations using
    output data

    :class:`REPOError` -- Raised by HDF5 repository interactions

    :class:`SYMMError` -- Raised by :mod:`opan.utils.symm` submodule functions

    :class:`VECTORError` -- Raised by :mod:`opan.utils.vector`
    submodule functions

    :class:`XYZError` -- Raised during parsing of or calculations using
    XYZ data


----------------

**API**

"""


# Module-level imports


class OPANError(Exception):
    """Base class for custom errors defined for OpenAnharmonic

    :class:`OPANError` is an abstract superclass of any custom errors built
    under the OpenAnharmonic umbrella. It defines all common methods shared
    among the various subtype error classes, such that the only contents that
    must be declared by a subclass are str class variables with contents
    identical to their names.  These are recognized by the
    :meth:`~OPANError.__metaclass__.__iter__`
    defined in :class:`~OPANError.__metaclass__` as being the
    set of valid typecodes.

    Arguments
    ---------
    tc  : str
        String representation of typecode to be associated with the
        :class:`OPANError` subclass instance. *Must* be a validly
        constructed typecode defined for the relevant subclass.
    msg : str
        Explanation of the nature of the error being reported.
    src : str
        Further detail of the code/file source of the error behavior.

    Raises
    ------
    NotImplementedError
        Upon attempt to instantiate abstract :class:`OPANError` base class.

    Attributes
    ----------
    msg : str
        Explanation of the nature of the error being reported.
    src : str
        Further detail of the code source of the error behavior.
    subclass_name   : str
        String representation of the name of the OPANError subclass of the
        instance.
    tc  : str
        String typecode associated with the instance


    .. class:: __metaclass__(type)

        Metaclass providing ability to iterate over typecodes.

        With this metaclass, iterating over the class itself (rather than an
        instance) yields the valid typecodes for the class.

        .. method:: __iter__()

            Iterate over all defined typecodes.

            Generator iterating over all class variables whose names match
            their contents. For a properly constructed
            :class:`~opan.error.OPANError`
            subclass, these are identical to the typecodes.

            **Example:**

            >>> 'xyzfile' in opan.error.XYZError
            True

    """

    def __init__(self, tc, msg, src):

        # Import(s)
        import re

        # Check for and complain at instantiation of base class.
        if type(self) == OPANError:
            raise(NotImplementedError("OPANError base class is abstract."))
        ## end if

        # Quick RegEx to extract the name of the subclass.
        self.subclass_name = re.search(self.__module__ + "\\.(?P<cls>\w+)'", \
                    str(self.__class__), re.I).group("cls")

        # Check for valid typecode and throw a more descriptive error if
        #  invalid.
        if not tc in self.__class__:
            raise(KeyError("Invalid " + self.subclass_name + \
                                                " typecode: " + str(tc)))
        ## end if

        # Store error content
        self.tc = tc
        self.msg = msg
        self.src = src

    ## end def __init__


    def __str__(self):  # pragma: no cover   (str rep has no code significance)
        """ String representation of an :class:`OPANError` subclass instance.

        Implemented primarily so that the error stack handling of the Python
        interpreter will provide useful information to the user.

        Return value is constructed as:

        ``(typecode string) Error message: Error source``

        Returns
        -------
        str
            String representation of the instance.

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

## end class OPANError


class XYZError(OPANError):
    """Error relating to parsing of or calculation using XYZ data.

    See the :class:`OPANError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Inconsistent geometry in an OpenBabel XYZ file
    #:
    #: * ORCA -- `.xyz` or `.trj`
    xyzfile = 'xyzfile'

    #: Object already initialized (overwrite not supported)
    overwrite = 'overwrite'

    #: Dihedral angle calculation requested for a set of atoms containing
    #: an insufficiently nonlinear trio of atoms
    dihed = 'dihed'

    #: Insufficient non-parallel character in some manner of calculation
    nonprl = 'nonprl'

## end class XYZError


class GRADError(OPANError):
    """Error relating to parsing of or calculation from gradient data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
            numats     : Invalid number of atoms specification,
                            or specification not found
            en         : Energy specification not found
            gradblock  : Invalid gradient block or gradient block
                            not found
            geomblock  : Invalid geometry block or geometry block not found

    """

    # Typecodes as class-level variables
    numats = 'numats'
    en = 'en'
    gradblock = 'gradblock'
    geomblock = 'geomblock'

## end class ENGRADError


class OUTPUTError(OPANError):
    """Error relating to parsing of or calculation from output data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
            {none yet}

    """

## end class OUTPUTError


class HESSError(OPANError):
    """Error relating to parsing of or calculation from Hessian data.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
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

    # Typecodes as class-level variables
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

## end class HESSError


class SYMMError(OPANError):
    """Error relating to determination of molecular symmetry.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
            notfound    : Symmetry element expected but not found.

    """

    # Typecodes as class-level variables
    notfound = 'notfound'

## end class SYMMError


class REPOError(OPANError):
    """Error relating to HDF5 repository interactions.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
            status      : HDF5 repo in improper status for requested operation.
            group       : Problem with a group in linked HDF5 file
            data        : Problem with a dataset in linked HDF5 file

    """

    # Typecodes as class-level variables
    status = 'status'
    data = 'data'
    group = 'group'

## end class REPOError


class ANHARMError(OPANError):
    """Error relating to :class:`~opan.anharm.OPAN_ANHARM` actions.

    See the OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
            repo        : OPAN_REPO conflict -- no repo bound when assignment
                            attempted, or attempt made to bind new repo when
                            one already bound
            status      : OPAN_ANHARM internal variables in inappropriate
                            status for the requested operation

    """

    # Typecodes as class-level variables
    repo = 'repo'
    status = 'status'

## end class ANHARMError


class INERTIAError(OPANError):
    """Error relating to utils.inertia submodule functions.

    See OPANError.__doc__ for more information.

    Attributes:
        tc, msg, src, subclass_name are inherited from OPANError

        typecodes:
            neg_moment  : A negative principal inertial moment was computed
            top_type    : No valid molecular top type was identified
            bad_geom    : A geometry being parsed was unsuitable for a
                            particular type of calculation/manipulation

    """

    # Typecodes as class-level variables
    neg_moment = 'neg_moment'
    top_type = 'top_type'
    bad_geom = 'bad_geom'

## end class INERTIAError


class VECTORError(OPANError):
    """Error relating to :mod:`opan.utils.vector` submodule functions.

    See the :class:`OPANError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Insufficient non-parallel character in some manner of calculation
    nonprl = 'nonprl'

    #: Vectors which should have been orthonormal were determined not to be
    orthonorm = 'orthonorm'

## end class VECTORError



if __name__ == '__main__':  # pragma: no cover
    print("Module not executable")
