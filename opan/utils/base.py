#-------------------------------------------------------------------------------
# Name:        utils
# Purpose:     Module containing utility functions for OpenAnharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     15 Aug 2014
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

"""General purpose utility functions for OpenAnharmonic.


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
    **NON-ITERABLE**, through the use of a local subclass of `str` that
    cripples the ``__iter__`` method. Any strings passed are returned
    in the packed tuples as standard, **ITERABLE** instances of `str`, however.

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
    tups :  `list` of `tuple`
        Number of tuples returned is equal to the length of the iterables
        passed in `*args`

    Raises
    ------
    ~exceptions.ValueError
        If any iterable objects are of different lengths

    """

    # Imports
    #import numpy as np

    # Debug flag
    _DEBUG = False

    # Uninitialized test value
    UNINIT_VAL = -1

    # Print the input if in debug mode
    if _DEBUG: # pragma: no cover
        print("args = " + str(args))

    # Non-iterable subclass of str
    class str_noiter(str):
        def __iter__(self):
            raise(NotImplementedError("Non-iterable string"))
        ## end def __iter__
    ## end class str_noiter

    # Re-wrap input arguments with non-iterable strings if required
    mod_args = range(len(args))
    for idx in range(len(args)):
        if isinstance(args[idx], str):
            mod_args[idx] = str_noiter(args[idx])
        else:
            mod_args[idx] = args[idx]
        ## end if
    ## next idx

    # Initialize the info variable for the iterable length and scan all
    #  of the modified args members for iterability and length
    iter_len = UNINIT_VAL
    for idx in range(len(mod_args)):
        if iterable(mod_args[idx]):
            if iter_len == UNINIT_VAL:
                iter_len = len(mod_args[idx])
            else:
                if not iter_len == len(mod_args[idx]):
                    raise(ValueError("All iterable items must be of " + \
                            "equal length."))
                ## end if
            ## end if
        ## end if
    ## next idx

    # If everything is non-iterable, just return the args tuple wrapped in
    #  a list
    if iter_len == UNINIT_VAL:
        return [args]
    ## end if

    # Some items are iterable, so must construct iteration strings
    # Initialize strings for the iteration variables, the zip contents,
    #  and the output block of the final tuple-building for-loop
    itstr = ""
    zipstr = ""
    callstr = ""

    # Append suitable blips as things are iterable. Iteration variables are
    #  of the form "a###", with ### as digit(s).  Variables being zipped are
    #  the various elements of ***args***, which if there are any strings
    #  present will be instances of the **base str class**, rather than the
    #  local non-iterable string class, which is the desired behavior.
    # The call string is built either by iteration over iterables, or by
    #  repeated placement of non-iterables.
    for idx in range(len(mod_args)):
        if iterable(mod_args[idx]):
            itstr = itstr + "a" + str(idx) + ", "
            zipstr = zipstr + "args[" + str(idx) + "], "
            callstr = callstr + "a" + str(idx) + ", "
        else:
            callstr = callstr + "args[" + str(idx) + "], "
        ## end if
    ## next idx

    # Looping is always required if code gets to this point
    # Trim the excess comma and space from the strings
    itstr = itstr[:len(itstr) - 2]
    zipstr = zipstr[:len(zipstr) - 2]
    callstr = callstr[:len(callstr) - 2]

    # Build the start of the evaluation string
    evalstr = "[ (" + callstr + ") for " + itstr + " in "

    # Only zip() if more than one thing is being looped over
    evalstr = evalstr + ("zip(" if itstr.count(",") >= 1 else "")

    # Provide the zip string
    evalstr = evalstr + zipstr

    # Only close the zip() if more than one thing is being looped
    evalstr = evalstr + (")]" if itstr.count(",") >= 1 else "]")

    # Dump the built strings if in debug mode
    if _DEBUG:  # pragma: no cover
        print("evalstr = " + evalstr)
        print("itstr = " + itstr)
        print("zipstr = " + zipstr)
        print("callstr = " + callstr)
    ## end if

    # eval() the for loop to obtain tuples of arguments
    tups = eval(evalstr)

    # Dump the resulting tuples, if in debug mode
    if _DEBUG:  # pragma: no cover
        print("tups = " + str(tups))
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
    int
        Value of Kronecker delta for provided indices, as tested by
        Python "\ `==`\ "

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
    invar   : arbitrary
        Value to be typecast.
    totype  : `type`
        Type to which `invar` is to be cast.

    Returns
    -------
    outvar  : `type 'totype'`
        Typecast version of `invar`

    Raises
    ------
    ~exceptions.TypeError
        If result of typecast is not of type `totype`

    """

    # Make the typecast. Just use Python built-in exceptioning
    outvar = totype(invar)

    # Check that the cast type matches
    if not isinstance(outvar, totype):
        raise(TypeError("Result of cast to '" + str(totype) + \
                    "' is type '" + str(type(outvar)) ))
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
    el_time : `int` or `float`
        Time interval in seconds to be converted to h/m/s format

    Returns
    -------
    str
        String timestamp in h/m/s format
    """

    # Calc hours
    hrs = el_time // 3600.0

    # Calc minutes
    mins = (el_time % 3600.0) // 60.0

    # Calc seconds
    secs = el_time % 60.0

    # Construct timestamp string
    stamp = str(int(hrs)) + 'h ' + str(int(mins)) + 'm ' + str(int(secs)) + 's'

    # Return
    return stamp

## end def make_timestamp


@_arraysqueeze(0,1,2,3)
def check_geom(c1, a1, c2, a2, tol=_DEF.XYZ_Coord_Match_Tol):
    """ Check for consistency of two geometries and atom symbol lists

    Cartesian coordinates are considered consistent with the input
    coords if each component matches to within `tol`.  If coords or
    atoms vectors are passed that are of mismatched lengths, a
    ``False`` value is returned.

    Both coords vectors must be three times the length of the atoms vectors
    or a :exc:`~exceptions.ValueError` is raised.

    Parameters
    ----------
    c1      : length-3N ``np.float_``
        Vector of first set of stacked 'lab-frame' Cartesian coordinates

    a1      : length-N `str` or `int`
        Vector of first set of atom symbols or atomic numbers

    c2      : length-3N ``np.float_``
        Vector of second set of stacked 'lab-frame' Cartesian coordinates

    a2      : length-N `str` or `int`
        Vector of second set of atom symbols or atomic numbers

    tol    : float, optional
        Tolerance for acceptable deviation of each geometry coordinate
        from that in the reference instance to still be considered
        matching. Default value is specified by
        :attr:`opan.const.DEF.XYZ_Coord_Match_Tol`)

    Returns
    -------
    match  : bool
        Whether input coords and atoms match (\ ``True``\ ) or
        not (\ ``False``\ )

    fail_type  : `str` or ``None``
        Type of check failure

        If `match` == ``True``:

            Returns ``None``

        If `match` == ``False``:

            A string code describing the reason for the failed match:

                `coord_dim_mismatch` -- Mismatch in coordinate vector sizes

                `atom_dim_mismatch`  -- Mismatch in atom vector sizes

                `coord_mismatch`     -- Mismatch in one or more coordinates

                `atom_mismatch`      -- Mismatch in one or more atoms

                **#TODO:** ``opan.utils.check_geom``: Convert ``fail_type`` to Enum

    fail_loc   : length-3N `bool` or length-N `bool` or ``None``
        Mismatched elements

        If `match` == ``True``:

            Returns ``None``

        If `match` == ``False``:

            For "array-level" problems such as a dimension mismatch, a
            ``None`` value is returned.

            For "element-level" problems, an ``np.array`` vector is returned
            indicating positions of mismatch in either `coords` or `atoms`,
            depending on the value of `fail_type`.

                ``True`` elements indicate **MATCHING** values

                ``False`` elements mark **MISMATCHES**

    Raises
    ------
    ValueError
        If array lengths are inconsistent ::

            if len(c1) != 3 * len(a1) or len(c2) != 3 * len(a2):
                raise(ValueError(...))

    """

    # Import(s)
    from ..const import atomNum
    import numpy as np

    # Initialize return value to success condition
    match = True

    #** Check coords for suitable shape. Assume 1-D np.arrays.
    if not len(c1.shape) == 1:
        # Cannot coerce to vector; complain.
        raise(ValueError(("'c1' is not a vector.")))
    ## end if
    if not len(c2.shape) == 1:
        # Cannot coerce to vector; complain.
        raise(ValueError(("'c2' is not a vector.")))
    ## end if

    #** Check atoms for suitable shape. Assume lists of strings, so
    # convert to np.array to check.
    if not len(a1.shape) == 1:
        # Not a vector; complain
        raise(ValueError(("'a1' is not a simple list.")))
    ## end if
    if not len(a2.shape) == 1:
        # Not a vector; complain.
        raise(ValueError(("'a2' is not a simple list.")))
    ## end if

    #** Confirm proper lengths of coords vs atoms
    if not c1.shape[0] == 3 * a1.shape[0]:
        raise(ValueError("len(c1) != 3*len(a1)"))
    ## end if
    if not c2.shape[0] == 3 * a2.shape[0]:
        raise(ValueError("len(c2) != 3*len(a2)"))
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
    #  instance atom arrays to atomNums to allow np.equals comparison.
    if np.issubdtype(a1.dtype, np.dtype('str')):
        # Presume atomic symbol data and attempt conversion
        a1 = np.array([atomNum[e] for e in a1])
    ## end if
    if np.issubdtype(a2.dtype, np.dtype('str')):
        # Presume atomic symbol data and attempt conversion
        a2 = np.array([atomNum[e] for e in a2])
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


def template_subst(template, subs, delims=['<', '>']):
    """ Perform substitution of content into tagged string.

    For substitutions into template input files for external computational
    packages, no checks for valid syntax are performed.

    The first element of each 2-tuple of `subs` corresponds to a delimited
    substitution tag to be replaced in `template` by the entire
    text of the second element. For example, the tuple ``("ABC", "text")`` would
    convert ``The <ABC> is working`` to  ``The text is working``, using the
    default delimiters of '<' and '>'. Substitutions are performed in
    iteration order from the `subs` iterable. Recursive substitution
    as the tag parsing proceeds is thus
    feasible if an ordered iterable is used.

    Start and end delimiters for the tags are modified by `delims`. For
    example, to substitute a tag of the form **{\|TAG\|}**, the tuple
    ``("{|","|}")`` should be passed to `subs_delims`.  Any elements in
    `delims` past the second are ignored. No checking is
    performed for whether the delimiters are "sensible" or not.

    Parameters
    ----------
    template : `str`
        Template containing tags delimited by `subs_delims`,
        with tag names and substitution contents provided in `subs`

    subs    : `tuple` of `2-tuples` of `str`
        Each `2-tuple` contains a tag name and corresponding content to be
        substituted into the provided template.

    delims : iterable of `str`
        Iterable containing the 'open' and 'close' strings used to mark tags
        in the template, which are drawn from elements zero and one,
        respectively. Any elements beyond these are ignored.

    Returns
    -------
    str
        String generated from the parsed template, with all tag
        substitutions performed.

    """

    # Store the template into the working variable
    input_text = template

    # Iterate over subs and perform the .replace() calls
    for tup in subs:
        input_text = input_text.replace(
                delims[0] + tup[0] + delims[1], tup[1])
    ## next tup

    # Return the result
    return input_text

## end def template_subst


def iterable(y):
    """Check whether or not an object supports iteration.

    Adapted directly from NumPy ~= 1.10 at commit `46d2e83
    <https://github.com/numpy/numpy/tree/
    46d2e8356760e7549d0c80da9fe232177924183c/numpy/lib/
    function_base.py#L48-L76>`__.

    Parameters
    ----------
    y : object
      Object to be tested.

    Returns
    -------
    bool
      Returns ``True`` if the object has an iterator method or is a sequence,
      and ``False`` otherwise.

    Examples
    --------
    >>> np.iterable([1, 2, 3])
    True
    >>> np.iterable(2)
    False

    """
    try:
        iter(y)
    except:
        return False
    return True


## end def iterable


if __name__ == '__main__': # pragma: no cover
    print("Module not executable.")





