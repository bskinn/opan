#-------------------------------------------------------------------------------
# Name:        utils
# Purpose:     Module containing utility functions for OpenAnharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     15 Aug 2014
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

"""General purpose utility functions for OpenAnharmonic.

*This module docstring is not used in the Sphinx docs.*

check_geom       -- Confirm two OpenBabel geometries (atom types and
                        coordinates) match to within a specified tolerance
delta_fxn        -- Generalized Kronecker delta function
make_timestamp   -- Construct a string time-elapsed timestamp in h/m/s format
pack_tups        -- Pack an arbitrary combination of iterables and non-
                        iterables into a list of tuples
safe_cast        -- Robustified casting with a post-check to confirm the cast
                        actually resulted in the proper type
template_subst   -- Perform a field-based substitution into a string template


"""

# Imports
from ..const import DEF as _DEF
from .decorate import arraysqueeze as _arraysqueeze


def pack_tups(*args):
    """Pack an arbitrary set of iterables and non-iterables into tuples.

    Function packs a set of inputs with arbitrary iterability into tuples.
    Iterability is tested with :func:`iterable`. Non-iterable inputs
    are repeated in each output tuple. Iterable inputs are expanded
    uniformly across the output tuples.  For consistency, all iterables must
    be the same length.

    The input arguments are parsed such that bare strings are treated as
    **NON-ITERABLE**, through the use of a local subclass of |str| that
    cripples the ``__iter__()`` method. Any strings passed are returned
    in the packed tuples as standard, **ITERABLE** instances of |str|, however.

    The order of the input arguments is retained within each output tuple.

    No structural conversion is attempted on the arguments.

    If all inputs are non-iterable, a list containing a single `tuple` will be
    returned.

    Parameters
    ----------
    \*args
        Arbitrary number of arbitrary mix of iterable and non-iterable
        objects to be packed into tuples.

    Returns
    -------
    tups
        |list| of |tuple| --
        Number of tuples returned is equal to the length of the iterables
        passed in ``*args``

    Raises
    ------
    ~exceptions.ValueError
        If any iterable objects are of different lengths

    """

    # Imports
    import numpy as np

    # Debug flag
    _DEBUG = False

    # Marker value for non-iterable items
    NOT_ITER = -1

    # Uninitialized test value
    UNINIT_VAL = -1

    # Print the input if in debug mode
    if _DEBUG: # pragma: no cover
        print("args = {0}".format(args))

    # Non-iterable subclass of str
    class StrNoIter(str):
        """ Non-iterable subclass of |str|. """
        def __iter__(self):
            raise NotImplementedError("Non-iterable string")
        ## end def __iter__
    ## end class StrNoIter

    # Re-wrap input arguments with non-iterable strings if required
    mod_args = [(StrNoIter(a) if isinstance(a, str) else a) for a in args]

    # Determine the length or non-iterable status of each item and store
    #  the maximum value (depends on NOT_ITER < 0)
    iterlens = [(len(a) if iterable(a) else NOT_ITER) for a in mod_args]
    maxiter = max(iterlens)

    # Check to ensure all iterables are the same length
    if not all(map(lambda v: v in (NOT_ITER, maxiter), iterlens)):
        raise ValueError("All iterable items must be of equal length")
    ## end if

    # If everything is non-iterable, just return the args tuple wrapped in
    #  a list (as above, depends on NOT_ITER < 0)
    if maxiter == NOT_ITER:
        return [args]
    ## end if

    # Swap any non-iterables for a suitable length repeat, and zip to
    #  tuples for return
    tups = zip(*[(np.repeat(a, maxiter) if l == NOT_ITER else a)
            for (a,l) in zip(mod_args, iterlens)])

    # Dump the resulting tuples, if in debug mode
    if _DEBUG:  # pragma: no cover
        print("tups = {0}".format(tups))
    ## end if

    # Return the tuples
    return tups

## end def pack_tups


def delta_fxn(a, b):
    """Kronecker delta for objects `a` and `b`.

    Parameters
    ----------
    a :
        First object
    b :
        Second object

    Returns
    -------
    delta
        |int| --
        Value of Kronecker delta for provided indices, as tested by
        Python "`==`"

    """

    return (1 if a == b else 0)

## end def delta_fxn


def safe_cast(invar, totype):
    """Performs a "safe" typecast.

    Ensures that `invar` properly casts to `totype`. Checks after
    casting that the result is actually of type `totype`. Any exceptions raised
    by the typecast itself are unhandled.

    Parameters
    ----------
    invar
        (arbitrary) -- Value to be typecast.

    totype
        |type| --  Type to which `invar` is to be cast.

    Returns
    -------
    outvar
        `type 'totype'` --  Typecast version of `invar`

    Raises
    ------
    ~exceptions.TypeError
        If result of typecast is not of type `totype`

    """

    # Make the typecast. Just use Python built-in exceptioning
    outvar = totype(invar)

    # Check that the cast type matches
    if not isinstance(outvar, totype):
        raise TypeError("Result of cast to '{0}' is '{1}'"
                                            .format(totype, type(outvar)))
    ## end if

    # Success; return the cast value
    return outvar

## end def safe_cast


def make_timestamp(el_time):
    """ Generate an hour-minutes-seconds timestamp from an interval in seconds.

    Assumes numeric input of a time interval in seconds.  Converts this
    interval to a string of the format "#h #m #s", indicating the number of
    hours, minutes, and seconds in the interval.  Intervals greater than 24h
    are unproblematic.

    Parameters
    ----------
    el_time
        |int| or |float| --
        Time interval in seconds to be converted to h/m/s format

    Returns
    -------
    stamp
        |str| -- String timestamp in #h #m #s format

    """

    # Calc hours
    hrs = el_time // 3600.0

    # Calc minutes
    mins = (el_time % 3600.0) // 60.0

    # Calc seconds
    secs = el_time % 60.0

    # Construct timestamp string
    stamp = "{0}h {1}m {2}s".format(int(hrs), int(mins), int(secs))

    # Return
    return stamp

## end def make_timestamp


@_arraysqueeze(0,1,2,3)
def check_geom(c1, a1, c2, a2, tol=_DEF.XYZ_COORD_MATCH_TOL):
    """ Check for consistency of two geometries and atom symbol lists

    Cartesian coordinates are considered consistent with the input
    coords if each component matches to within `tol`.  If coords or
    atoms vectors are passed that are of mismatched lengths, a
    |False| value is returned.

    Both coords vectors must be three times the length of the atoms vectors
    or a :exc:`~exceptions.ValueError` is raised.

    Parameters
    ----------
    c1
        length-3N |npfloat|_ --
        Vector of first set of stacked 'lab-frame' Cartesian coordinates

    a1
        length-N |str| or |int| --
        Vector of first set of atom symbols or atomic numbers

    c2
        length-3N |npfloat|_ --
        Vector of second set of stacked 'lab-frame' Cartesian coordinates

    a2
        length-N |str| or |int| --
        Vector of second set of atom symbols or atomic numbers

    tol
        |float|, optional --
        Tolerance for acceptable deviation of each geometry coordinate
        from that in the reference instance to still be considered
        matching. Default value is specified by
        :attr:`opan.const.DEF.XYZ_COORD_MATCH_TOL`)

    Returns
    -------
    match
        |bool| --
        Whether input coords and atoms match (|True|) or
        not (|False|)

    fail_type
        |str| or |None| -- Type of check failure

        If `match` == |True|:

            Returns as |None|

        If `match` == |False|:

            A string code describing the reason for the failed match:

                `coord_dim_mismatch` -- Mismatch in coordinate vector sizes

                `atom_dim_mismatch`  -- Mismatch in atom vector sizes

                `coord_mismatch`     -- Mismatch in one or more coordinates

                `atom_mismatch`      -- Mismatch in one or more atoms

                **#TODO:** :func:`~opan.utils.base.check_geom`: Convert `fail_type` to Enum

    fail_loc
        length-3N |bool| or length-N |bool| or |None| --
        Mismatched elements

        If `match` == |True|:

            Returns as |None|

        If `match` == |False|:

            For "array-level" problems such as a dimension mismatch, a
            |None| value is returned.

            For "element-level" problems, a vector is returned
            indicating positions of mismatch in either `coords` or `atoms`,
            depending on the value of `fail_type`.

                |True| elements indicate **MATCHING** values

                |False| elements mark **MISMATCHES**

    Raises
    ------
    ~exceptions.ValueError
        If a pair of coords & atoms array lengths is inconsistent ::

            if len(c1) != 3 * len(a1) or len(c2) != 3 * len(a2):
                raise (ValueError(...)

    """

    # Import(s)
    from ..const import atom_num
    import numpy as np

    # Initialize return value to success condition
    match = True

    #** Check coords for suitable shape. Assume 1-D np.arrays.
    if not len(c1.shape) == 1:
        # Cannot coerce to vector; complain.
        raise ValueError(("'c1' is not a vector."))
    ## end if
    if not len(c2.shape) == 1:
        # Cannot coerce to vector; complain.
        raise ValueError(("'c2' is not a vector."))
    ## end if

    #** Check atoms for suitable shape. Assume lists of strings, so
    # convert to np.array to check.
    if not len(a1.shape) == 1:
        # Not a vector; complain
        raise ValueError(("'a1' is not a simple list."))
    ## end if
    if not len(a2.shape) == 1:
        # Not a vector; complain.
        raise ValueError(("'a2' is not a simple list."))
    ## end if

    #** Confirm proper lengths of coords vs atoms
    if not c1.shape[0] == 3 * a1.shape[0]:
        raise ValueError("len(c1) != 3*len(a1)")
    ## end if
    if not c2.shape[0] == 3 * a2.shape[0]:
        raise ValueError("len(c2) != 3*len(a2)")
    ## end if

    #** Confirm matching lengths of coords and atoms w/corresponding
    #  objects among the two geometries
    if not c1.shape[0] == c2.shape[0]:
        match = False
        fail_type = "coord_dim_mismatch"
        return match, fail_type, None
    ## end if
    if not a1.shape[0] == a2.shape[0]:
        match = False
        fail_type = "atom_dim_mismatch"
        return match, fail_type, None
    ## end if

    #** Element-wise check for geometry match to within 'tol'
    fail_loc = np.less_equal(np.abs(np.subtract(c1,c2)), tol)
    if sum(fail_loc) != c2.shape[0]:
        # Count of matching coordinates should equal the number of
        #  coordinates. If not, complain with 'coord_mismatch' fail type.
        match = False
        fail_type = "coord_mismatch"
        return match, fail_type, fail_loc
    ## end if

    #** Element-wise check for atoms match. Quietly convert both input and
    #  instance atom arrays to atom_nums to allow np.equals comparison.
    if np.issubdtype(a1.dtype, np.dtype('str')):
        # Presume atomic symbol data and attempt conversion
        a1 = np.array([atom_num[e] for e in a1])
    ## end if
    if np.issubdtype(a2.dtype, np.dtype('str')):
        # Presume atomic symbol data and attempt conversion
        a2 = np.array([atom_num[e] for e in a2])
    ## end if
    fail_loc = np.equal(a1, a2)

    #** Perform the test to ensure all atoms match.
    if sum(fail_loc) != a2.shape[0]:
        # Count of matching atoms should equal number of atoms. If not,
        #  complain with the 'atom_mismatch' fail type.
        match = False
        fail_type = "atom_mismatch"
        return match, fail_type, fail_loc

    #** If reached here, all tests passed; return success.
    return match, None, None

## end def check_geom


def template_subst(template, subs, delims=('<', '>')):
    """ Perform substitution of content into tagged string.

    For substitutions into template input files for external computational
    packages, no checks for valid syntax are performed.

    Each key in `subs` corresponds to a delimited
    substitution tag to be replaced in `template` by the entire text of the
    value of that key. For example, the dict ``{"ABC": "text"}`` would
    convert ``The <ABC> is working`` to  ``The text is working``, using the
    default delimiters of '<' and '>'. Substitutions are performed in
    iteration order from `subs`; recursive substitution
    as the tag parsing proceeds is thus
    feasible if an :class:`~collections.OrderedDict` is used and substitution
    key/value pairs are added in the proper order.

    Start and end delimiters for the tags are modified by `delims`. For
    example, to substitute a tag of the form **{\|TAG\|}**, the tuple
    ``("{|","|}")`` should be passed to `subs_delims`.  Any elements in
    `delims` past the second are ignored. No checking is
    performed for whether the delimiters are "sensible" or not.

    Parameters
    ----------
    template
        |str| --
        Template containing tags delimited by `subs_delims`,
        with tag names and substitution contents provided in `subs`

    subs
        |dict| of |str| --
        Each item's key and value are the tag name and corresponding content to
        be substituted into the provided template.

    delims
        iterable of |str| --
        Iterable containing the 'open' and 'close' strings used to mark tags
        in the template, which are drawn from elements zero and one,
        respectively. Any elements beyond these are ignored.

    Returns
    -------
    subst_text
        |str| --
        String generated from the parsed template, with all tag
        substitutions performed.

    """

    # Store the template into the working variable
    subst_text = template

    # Iterate over subs and perform the .replace() calls
    for (k,v) in subs:
        subst_text = subst_text.replace(
                delims[0] + k + delims[1], v)
    ## next tup

    # Return the result
    return subst_text

## end def template_subst


def iterable(y):
    """Check whether or not an object supports iteration.

    Adapted directly from NumPy ~= 1.10 at commit `46d2e83
    <https://github.com/numpy/numpy/tree/
    46d2e8356760e7549d0c80da9fe232177924183c/numpy/lib/
    function_base.py#L48-L76>`__.

    Parameters
    ----------
    y
        (arbitrary) -- Object to be tested.

    Returns
    -------
    test
        |bool| --
        Returns |False| if :func:`iter` raises an exception when `y` is
        passed to it; |True| otherwise.

    Examples
    --------
    >>> np.iterable([1, 2, 3])
    True
    >>> np.iterable(2)
    False

    """
    try:
        iter(y)
    except Exception:
        return False
    return True


## end def iterable


def assert_npfloatarray(obj, varname, desc, exc, tc, errsrc):
    """ Assert a value is an |nparray| of NumPy floats.

    Pass |None| to `varname` if `obj` itself is to be checked.
    Otherwise, `varname` is the string name of the attribute of `obj` to
    check.  In either case, `desc` is a string description of the
    object to be checked, for use in raising of exceptions.

    Raises the exception `exc` with typecode `tc` if the indicated
    object is determined not to be an |nparray|, with a NumPy float dtype.

    Intended primarily to serve as an early check for
    proper implementation of subclasses of
    :class:`~opan.grad.SuperOpanGrad` and
    :class:`~opan.hess.SuperOpanHess`. Early type-checking of key
    attributes will hopefully avoid confusing bugs downstream.

    Parameters
    ----------
    obj
        (arbitrary) --
        Object to be checked, or object with attribute to be checked.

    varname
        |str| or |None| --
        Name of the attribute of `obj` to be type-checked. |None|
        indicates to check `obj` itself.

    desc
        |str| --
        Description of the object being checked to be used in any
        raised exceptions.

    exc
        Subclass of :class:`~opan.error.OpanError` to be raised on
        a failed typecheck.

    tc
        Typecode of `exc` to be raised on a failed typecheck.

    errsrc
        |str| --
        String description of the source of the data leading to a
        failed typecheck.

    """

    # Imports
    import numpy as np

    # Check for whether member or object is to be checked
    if varname is None:
        var = obj
    else:
        # Try to get the variable to be typechecked
        try:
            var = getattr(obj, varname)
        except AttributeError:
            raise exc(tc, "Attribute '{0}' not defined in '{1}'"
                    .format(varname, obj), errsrc)
        ## end try
    ## end if

    # Try to pull the np dtype off of it
    try:
        dt = var.dtype
    except AttributeError:
        raise exc(tc, "'{0}' is not an np.array (lacks a 'dtype' member)"
                    .format(desc), errsrc)
    ## end try

    # Confirm dtype inherits from np.float
    if not np.issubdtype(dt, np.float):
        raise exc(tc, "'{0}' is not an np.array of np.float".format(desc),
                errsrc)
    ## end if

## end def assert_npfloatarray



if __name__ == '__main__': # pragma: no cover
    print("Module not executable.")





