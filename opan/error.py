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

Error classes are subclassed from :exc:`~exceptions.Exception` via an
abstract superclass, :class:`~opan.error.OpanError`, which defines
several common features:

 * Storage of a 'typecode' and a 'source' string along with the error message
   to allow passing of more fine-grained information to the exception stack
 * Implementation of a :attr:`~opan.error.OpanError.__metaclass__` on
   :class:`OpanError` enabling typecode validity checking with ``is``
 * Re-implementation of :func:`~OpanError.__str__` to enhance the usefulness
   of stack messages when one of these errors is raised

:class:`OpanError` **Subclasses**

    :class:`AnharmError` -- Raised as a result of
    :class:`~opan.anharm.base.OpanAnharm` actions

    :class:`GradError` -- Raised during parsing of or calculations using
    gradient data

    :class:`HessError` -- Raised during parsing of or calculations using
    Hessian data

    :class:`InertiaError` -- Raised by :mod:`opan.utils.inertia`
    submodule functions

    :class:`OutputError` -- Raised during parsing of or calculations using
    output data

    :class:`RepoError` -- Raised by HDF5 repository interactions

    :class:`SymmError` -- Raised by :mod:`opan.utils.symm` submodule functions

    :class:`VectorError` -- Raised by :mod:`opan.utils.vector`
    submodule functions

    :class:`XYZError` -- Raised during parsing of or calculations using
    XYZ data


----------------

**API**

"""


# Module-level imports


class OpanError(Exception):
    """Base class for custom errors defined for OpenAnharmonic

    :class:`OpanError` is an abstract superclass of any custom errors defined
    under the OpenAnharmonic umbrella. It defines all common methods shared
    among the various subtype error classes, such that the only contents that
    must be declared by a subclass are |str| class variables with contents
    identical to their names.  These are recognized by the
    :meth:`~OpanError.__metaclass__.__iter__`
    defined in :class:`~OpanError.__metaclass__` as being the
    set of valid typecodes.

    Arguments
    ---------
    tc
        |str| -- String representation of typecode to be associated with the
        :class:`OpanError` subclass instance. *Must* be a validly
        constructed typecode defined for the relevant subclass.

    msg
        |str| -- Explanation of the nature of the error being reported

    src
        |str| -- Further detail of the code/file source of the error
        behavior, if relevant

    Raises
    ------
    ~exceptions.NotImplementedError
        Upon attempt to instantiate abstract :class:`OpanError` base class

    ~exceptions.KeyError
        Upon instantiation with an invalid typecode

    Attributes
    ----------
    msg
        |str| -- Explanation of the nature of the error being reported

    src
        |str| -- Further detail of the code source of the error behavior

    subclass_name
        |str| -- String representation of the :class:`OpanError` subclass name

    tc
        |str| -- String typecode associated with the instance


    .. class:: __metaclass__(type)

        Metaclass providing ability to iterate over typecodes.

        With this metaclass, iterating over the class itself (rather than an
        instance) yields the valid typecodes for the class.

        .. method:: __iter__()

            Iterate over all defined typecodes.

            Generator iterating over all class variables whose names match
            their contents. For a properly constructed
            :class:`~opan.error.OpanError`
            subclass, these are identical to the typecodes.

            **Example:**

            >>> 'XYZFILE' in opan.error.XYZError
            True

    """

    def __init__(self, tc, msg, src):

        # Import(s)
        import re

        # Check for and complain at instantiation of base class.
        if type(self) == OpanError:
            raise(NotImplementedError("OpanError base class is abstract."))
        ## end if

        # Quick RegEx to extract the name of the subclass.
        self.subclass_name = re.search(self.__module__ + "\\.(?P<cls>\w+)'",
                    str(self.__class__), re.I).group("cls")

        # Check for valid typecode and throw a more descriptive error if
        #  invalid.
        if not tc in self.__class__:
            raise(KeyError("Invalid {0} typecode: {1}"
                                            .format(self.subclass_name, tc)))
        ## end if

        # Store error content
        self.tc = tc
        self.msg = msg
        self.src = src

    ## end def __init__


    def __str__(self):  # pragma: no cover   (str rep has no code significance)
        """ String representation of an :class:`OpanError` subclass instance.

        Implemented primarily so that the error stack handling of the Python
        interpreter will provide useful information to the user.

        Return value is constructed as:

        ``(typecode string) Error message: Error source``

        Returns
        -------
        retstr
            |str| -- String representation of the instance.

        """

        # Store and return the descriptive string
        retstr = "({0}) {1}".format(self.tc, self.msg) + (
                            ": {0}".format(self.src) if self.src else "")
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

## end class OpanError


class XYZError(OpanError):
    """Error relating to parsing of or calculation using XYZ data.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Inconsistent geometry in an OpenBabel XYZ file
    #:
    #: * |orca| -- `.xyz` or `.trj`
    XYZFILE = 'XYZFILE'

    #: Object already initialized (overwrite not supported)
    OVERWRITE = 'OVERWRITE'

    #: Dihedral angle calculation requested for a set of atoms containing
    #: an insufficiently nonlinear trio of atoms
    DIHED = 'DIHED'

    #: Insufficient non-parallel character in some manner of calculation
    NONPRL = 'NONPRL'

## end class XYZError


class GradError(OpanError):
    """Error relating to parsing of or calculation from gradient data.

    Not all typecodes may be relevant for all software packages.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Invalid number-of-atoms specification, or specification not found
    NUMATS = 'NUMATS'

    #: Energy value not found
    ENERGY = 'ENERGY'

    #: Malformed or missing gradient block
    GRADBLOCK = 'GRADBLOCK'

    #: Malformed or missing geometry block
    GEOMBLOCK = 'GEOMBLOCK'

    #: Missing or invalid gradient data; SHOULD only be used by
    #: :class:`~opan.grad.SuperOpanGrad`
    BADGRAD = 'BADGRAD'

    #: Missing or invalid geometry data; SHOULD only be used by
    #: :class:`~opan.grad.SuperOpanGrad`
    BADGEOM = 'BADGEOM'

    #: Missing or invalid atom symbols; SHOULD only be used by
    #: :class:`~opan.grad.SuperOpanGrad`
    BADATOM = 'BADATOM'

    #: Object already initialized (overwrite not supported)
    OVERWRITE = 'OVERWRITE'

## end class GradError


class OutputError(OpanError):
    """Error relating to parsing of or calculation from output data.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    *(none yet)*

    """

## end class OutputError


class HessError(OpanError):
    """Error relating to parsing of or calculation from Hessian data.

    Not all typecodes may be relevant for all software packages.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Malformed or missing atom/geometry specification block
    AT_BLOCK = "AT_BLOCK"

    #: Malformed or missing Hessian block
    HESS_BLOCK = "HESS_BLOCK"

    #: Malformed or missing frequencies block
    FREQ_BLOCK = "FREQ_BLOCK"

    #: Malformed or missing normal modes block
    MODES_BLOCK = "MODES_BLOCK"

    #: Malformed dipole derivatives block
    DIPDER_BLOCK = "DIPDER_BLOCK"

    #: Malformed IR spectrum block
    IR_BLOCK = "IR_BLOCK"

    #: Malformed polarizability derivatives block
    POLDER_BLOCK = "POLDER_BLOCK"

    #: Malformed Raman spectrum block
    RAMAN_BLOCK = "RAMAN_BLOCK"

    #: Malformed job list block
    JOB_BLOCK = "JOB_BLOCK"

    #: Malformed mass-weighted-Hessian eigenvalues block
    EIGVAL_BLOCK = "EIGVAL_BLOCK"

    #: Malformed mass-weighted-Hessian eigenvectors block
    EIGVEC_BLOCK = "EIGVEC_BLOCK"

    #: Malformed or missing energy value
    ENERGY = "ENERGY"

    #: Malformed or missing temperature value
    TEMP = "TEMP"

    #: Missing or invalid Hessian data; SHOULD only be used by
    #: :class:`~opan.hess.SuperOpanHess`
    BADHESS = 'BADHESS'

    #: Missing or invalid geometry data; SHOULD only be used by
    #: :class:`~opan.hess.SuperOpanHess`
    BADGEOM = 'BADGEOM'

    #: Missing or invalid atom symbols; SHOULD only be used by
    #: :class:`~opan.hess.SuperOpanHess`
    BADATOM = 'BADATOM'

    #: Object already initialized (overwrite not supported)
    OVERWRITE = 'OVERWRITE'

## end class HessError


class SymmError(OpanError):
    """Error relating to :mod:`opan.utils.symm` submodule functions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #DOC: Add note about this being subject to future development?

    #: Symmetry element expected but not found
    NOTFOUND = 'NOTFOUND'

## end class SymmError


class RepoError(OpanError):
    """Error relating to HDF5 repository interactions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: HDF5 repo in improper status for requested operation
    STATUS = 'STATUS'

    #: Problem with a dataset in linked HDF5 :class:`h5py:File`
    DATA = 'DATA'

    #: Problem with a group in linked HDF5 :class:`h5py:File`
    GROUP = 'GROUP'

## end class RepoError


class AnharmError(OpanError):
    """Error relating to :class:`~opan.anharm.base.OpanAnharm` actions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #DOC: Add object references once the DOM is more established?

    #: OpanAnharmRepo conflict -- no repo bound when assignment attempted,
    #: or attempt made to bind new repo when one already bound
    REPO = 'REPO'

    #: :class:`~opan.anharm.base.OpanAnharm` internal variables in
    #: inappropriate state for the requested operation
    STATUS = 'STATUS'

## end class AnharmError


class InertiaError(OpanError):
    """Error relating to :mod:`opan.utils.inertia` submodule functions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: A negative principal inertial moment was computed
    NEG_MOMENT = 'NEG_MOMENT'

    #: No valid molecular top type was identified
    TOP_TYPE = 'TOP_TYPE'

    #: A geometry being parsed was unsuitable for a particular type of
    #: calculation/manipulation
    BAD_GEOM = 'BAD_GEOM'

## end class InertiaError


class VectorError(OpanError):
    """Error relating to :mod:`opan.utils.vector` submodule functions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Insufficient non-parallel character in some manner of calculation
    NONPRL = 'NONPRL'

    #: Vectors which should have been orthonormal were determined not to be
    ORTHONORM = 'ORTHONORM'

## end class VectorError



if __name__ == '__main__':  # pragma: no cover
    print("Module not executable")
