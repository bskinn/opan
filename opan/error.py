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
   :class:`OpanError` enabling typecode validity checking with `is`
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
    must be declared by a subclass are `str` class variables with contents
    identical to their names.  These are recognized by the
    :meth:`~OpanError.__metaclass__.__iter__`
    defined in :class:`~OpanError.__metaclass__` as being the
    set of valid typecodes.

    Arguments
    ---------
    tc  : str
        String representation of typecode to be associated with the
        :class:`OpanError` subclass instance. *Must* be a validly
        constructed typecode defined for the relevant subclass.
    msg : str
        Explanation of the nature of the error being reported
    src : str
        Further detail of the code/file source of the error behavior,
        if relevant

    Raises
    ------
    ~exceptions.NotImplementedError
        Upon attempt to instantiate abstract :class:`OpanError` base class.
    ~exceptions.KeyError
        Upon instantiation with an invalid typecode

    Attributes
    ----------
    msg : str
        Explanation of the nature of the error being reported
    src : str
        Further detail of the code source of the error behavior
    subclass_name   : str
        String representation of the :class:`OpanError` subclass name
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
            :class:`~opan.error.OpanError`
            subclass, these are identical to the typecodes.

            **Example:**

            >>> 'xyzfile' in opan.error.XYZError
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
        str
            String representation of the instance.

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
    xyzfile = 'xyzfile'

    #: Object already initialized (overwrite not supported)
    overwrite = 'overwrite'

    #: Dihedral angle calculation requested for a set of atoms containing
    #: an insufficiently nonlinear trio of atoms
    dihed = 'dihed'

    #: Insufficient non-parallel character in some manner of calculation
    nonprl = 'nonprl'

## end class XYZError


class GradError(OpanError):
    """Error relating to parsing of or calculation from gradient data.

    Not all typecodes may be relevant for all software packages.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Invalid number-of-atoms specification, or specification not found
    numats = 'numats'

    #: Energy value not found
    en = 'en'

    #: Malformed or missing gradient block
    gradblock = 'gradblock'

    #: Malformed or missing geometry block
    geomblock = 'geomblock'

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
    at_block = "at_block"

    #: Malformed or missing Hessian block
    hess_block = "hess_block"

    #: Malformed or missing frequencies block
    freq_block = "freq_block"

    #: Malformed or missing normal modes block
    modes_block = "modes_block"

    #: Malformed dipole derivatives block
    dipder_block = "dipder_block"

    #: Malformed IR spectrum block
    ir_block = "ir_block"

    #: Malformed polarizability derivatives block
    polder_block = "polder_block"

    #: Malformed Raman spectrum block
    raman_block = "raman_block"

    #: Malformed job list block
    job_block = "job_block"

    #: Malformed mass-weighted-Hessian eigenvalues block
    eigval_block = "eigval_block"

    #: Malformed mass-weighted-Hessian eigenvectors block
    eigvec_block = "eigvec_block"

    #: Malformed or missing energy value
    energy = "energy"

    #: Malformed or missing temperature value
    temp = "temp"

## end class HessError


class SymmError(OpanError):
    """Error relating to :mod:`opan.utils.symm` submodule functions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Symmetry element expected but not found
    notfound = 'notfound'

## end class SymmError


class RepoError(OpanError):
    """Error relating to HDF5 repository interactions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: HDF5 repo in improper status for requested operation
    status = 'status'

    #: Problem with a dataset in linked HDF5 :class:`h5py:File`
    data = 'data'

    #: Problem with a group in linked HDF5 :class:`h5py:File`
    group = 'group'

## end class RepoError


class AnharmError(OpanError):
    """Error relating to :class:`~opan.anharm.base.OpanAnharm` actions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: OpanAnharmRepo conflict -- no repo bound when assignment attempted,
    #: or attempt made to bind new repo when one already bound
    repo = 'repo'

    #: :class:`~opan.anharm.base.OpanAnharm` internal variables in
    #: inappropriate state for the requested operation
    status = 'status'

## end class AnharmError


class InertiaError(OpanError):
    """Error relating to :mod:`opan.utils.inertia` submodule functions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: A negative principal inertial moment was computed
    neg_moment = 'neg_moment'

    #: No valid molecular top type was identified
    top_type = 'top_type'

    #: A geometry being parsed was unsuitable for a particular type of
    #: calculation/manipulation
    bad_geom = 'bad_geom'

## end class InertiaError


class VectorError(OpanError):
    """Error relating to :mod:`opan.utils.vector` submodule functions.

    See the :class:`OpanError` documentation for more information on
    attributes, methods, etc.

    **Typecodes**

    """

    #: Insufficient non-parallel character in some manner of calculation
    nonprl = 'nonprl'

    #: Vectors which should have been orthonormal were determined not to be
    orthonorm = 'orthonorm'

## end class VectorError



if __name__ == '__main__':  # pragma: no cover
    print("Module not executable")
