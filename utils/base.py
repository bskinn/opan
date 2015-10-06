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


def pack_tups(*args):
    """Pack an arbitrary set of iterables and non-iterables into tuples.

    Function packs a set of inputs with arbitrary iterability into tuples.
        Non-iterable inputs are repeated in each output tuple. Iterable inputs
        are expanded uniformly across the output tuples.  For consistency, all
        iterables must be the same length.

    The order of the input arguments is retained within each output tuple.

    No structural conversion is attempted on the individual args.

    If all inputs are non-iterable, a list containing a single tuple will be
        returned.

    Parameters
    ----------
    args : just about anything?
        Arbitrary number of arbitrary mix of iterable and non-iterable
        objects to be packed into tuples.
        DOES NOT WORK RELIABLY on bare strings passed in
        #!TODO: pack_tups: Robustify to work properly on bare strings.

    Returns
    -------
    tups :  list of tuples
        Number of tuples returned is equal to the length of the iterables
        passed in args

    Raises
    ------
    ValueError :  If all iterable objects are not the same length
    """

    # Imports
    import numpy as np

    # Debug flag
    DEBUG = False

    # Uninitialized test value
    UNINIT_VAL = -1

    # Print the input if in debug mode
    if DEBUG:
        print("args = " + str(args))

    # Initialize the info variable for the iterable length and scan all
    #  args members for iterability and length
    iter_len = UNINIT_VAL
    for idx in range(len(args)):
        if np.iterable(args[idx]):
            if iter_len == UNINIT_VAL:
                iter_len = len(args[idx])
            else:
                if not iter_len == len(args[idx]):
                    raise(ValueError("All iterable items must be of " + \
                            "equal length."))
                ## end if
            ## end if
        ## end if
    ## next idx

    # If everything is non-iterable, just return the args tuple
    if iter_len == UNINIT_VAL:
        return args
    ## end if

    # Some items are iterable, so must construct iteration strings
    # Initialize strings for the iteration variables, the zip contents,
    #  and the output block of the final tuple-building for-loop
    itstr = ""
    zipstr = ""
    callstr = ""

    # Append suitable blips as things are iterable. Iteration variables are
    #  of the form "a###", with ### as digit(s).  Variables being zipped are
    #  the various elements of args.  The call string is built either by
    #  iteration over iterables, or by repeated placement of non-iterables.
    for idx in range(len(args)):
        if np.iterable(args[idx]):
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
    if DEBUG:
        print("evalstr = " + evalstr)
        print("itstr = " + itstr)
        print("zipstr = " + zipstr)
        print("callstr = " + callstr)
    ## end if

    # eval() the for loop to obtain tuples of arguments
    tups = eval(evalstr)

    # Dump the resulting tuples, if in debug mode
    if DEBUG:
        print("tups = " + str(tups))
    ## end if

    # Return the tuples
    return tups

## end def pack_tups


def delta_fxn(a, b):
    """Returns Kronecker delta for indices a and b.

    Parameters
    ----------
    a : scalar
        First index
    b : scalar
        Second index

    Returns
    -------
    delta : int
        Value of Kronecker delta for provided indices

    Raises
    ------
    TypeError : if nonscalar inputs are provided

    """

    # Imports
    from numpy import isscalar

    # Throw exception if either of a or b is not a scalar numeric
    # isscalar() apparently checks for non-matrix AND numeric...?
    if not isscalar(a):
        raise TypeError("orca_utils.delta_fxn() requires scalar "
                "numeric inputs" + chr(10) + "            "
                "Type of 'a' is: " + str(type(a)))

    if not isscalar(b):
        raise TypeError("orca_utils.delta_fxn() requires scalar "
                "numeric inputs" + chr(10) + "            "
                "Type of 'b' is: " + str(type(b)))

    return (1 if a == b else 0)

## end def delta_fxn


def safe_cast(invar, totype, check=True):
    """Performs a safe typecast.

    Ensures that 'invar' is castable to 'totype'. Optionally 'check' after
    casting that the result is actually of type 'totype'.

    Parameters
    ----------
    invar   : arbitrary
        Value to be typecast.
    totype  : type or method performing typecast
        Type to which 'invar' is to be cast.
    check   : bool
        If True, confirm that result from typecast is actually type 'totype'

    Returns
    -------
    outvar  : type 'totype'
        Typecast version of 'invar'

    Raises
    ------
    TypeError : If 'invar' does not typecast to 'totype'.  Also, if check==True,
        if result of typecast is not of type 'totype'.
    TypeError : If 'totype' is anything other than type 'type'

    """

    # Confirm 'totype' is a type
    if not isinstance(totype, type):
        raise(TypeError("'totype' is not a type object."))
    ## end if

    # Attempt the typecast. Just use Python built-in exceptioning
    outvar = totype(invar)

    # If indicated, check that the cast type matches
    if check:
        if not isinstance(outvar, totype):
            raise(TypeError("Result of cast to '" + str(totype) + \
                        "' is type '" + str(type(outvar)) + ".'"))
        ## end if
    ## end if

    # Success; return the cast value
    return outvar


def make_timestamp(el_time):
    """ Generate an hour-minutes-seconds timestamp from an interval in seconds.

    Assumes numeric input of a time interval in seconds.  Converts this
    interval to a string of the format "#h #m #s", indicating the number of
    hours, minutes, and seconds in the interval.  Intervals greater than 24h
    are unproblematic.

    Parameters
    ----------
    el_time : int or float
        Time interval to be converted to h/m/s format

    Returns
    -------
    stamp   : str
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


def check_geom(c1, a1, c2, a2, tol=_DEF.XYZ_Coord_Match_Tol):
    """ Check for consistency of two geometries and atom symbol lists

    Cartesian coordinates are considered consistent with the input
        coords if each component matches to within 'tol' (default value
        specified by orca_const._DEF.Coord_Match_Tol).  If coords or
        atoms vectors are passed that are of mismatched lengths, a
        False value is returned.

    Both coords vectors must be three times the length of the atoms vectors
        or a ValueError is raised.

    Parameters
    ----------
    c1      : 3N x 1 float
        Vector of first set of stacked 'lab-frame' Cartesian coordinates
    a1      : N x 1 string or int
        Vector of first set of atom symbols or atomic numbers
    c2      : 3N x 1 float
        Vector of second set of stacked 'lab-frame' Cartesian coordinates
    a2      : N x 1 string or int
        Vector of second set of atom symbols or atomic numbers
    tol    : float, optional
        Tolerance for acceptable deviation of each geometry coordinate
        from that in the reference instance to still be considered
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
            #TODO: orca_utils.check_geom: Convert fail_type to an Enum
    fail_loc   : 3N x 1 bool or N x 1 bool
        np.matrix() column vector indicating positions of mismatch in
        either coords or atoms, depending on the value of fail_type.
        True elements indicate corresponding *MATCHING* values; False
        elements mark *MISMATCHES*.

    Raises
    ------
    ValueError : If len(c#) != 3 * len(a#)
    """

    # Import(s)
    from ..const import atomNum
    import numpy as np

    # Initialize return value to success condition
    match = True

    #** Check coords for suitable shape. Must convert to np.array() since
    #  a np.matrix() vector will not reduce to a single dimension.
    if not len(np.array(c1).squeeze().shape) == 1:
        # Cannot coerce to vector; complain.
        raise(ValueError(("'c1' is not a vector.")))
    else:
        # Coercible to a vector. Store a copy as an np.matrix() vector
        in_g1 = np.matrix(np.array(c1).squeeze().copy()).transpose()
    ## end if

    if not len(np.array(c2).squeeze().shape) == 1:
        # Cannot coerce to vector; complain.
        raise(ValueError(("'c2' is not a vector.")))
    else:
        # Coercible to a vector. Store a copy as an np.matrix() vector
        in_g2 = np.matrix(np.array(c2).squeeze().copy()).transpose()
    ## end if

    #** Check atoms for suitable shape. Must convert to np.array() since
    #  a np.matrix() vector will not reduce to a single dimension.
    if not len(np.array(a1).squeeze().shape) == 1:
        # Cannot coerce to vector; complain.
        raise(ValueError(("'a1' is not a vector.")))
    else:
        # Coercible to a vector. Store a copy as an np.matrix() vector
        in_a1 = np.matrix(np.array(a1).squeeze().copy()).transpose()
    ## end if

    if not len(np.array(a2).squeeze().shape) == 1:
        # Cannot coerce to vector; complain.
        raise(ValueError(("'a2' is not a vector.")))
    else:
        # Coercible to a vector. Store a copy as an np.matrix() vector
        in_a2 = np.matrix(np.array(a2).squeeze().copy()).transpose()
    ## end if

    #** Confirm proper lengths of coords vs atoms
    if not in_g1.shape[0] == 3 * in_a1.shape[0]:
        raise(ValueError("len(c1) != 3*len(a1)"))
    ## end if
    if not in_g2.shape[0] == 3 * in_a2.shape[0]:
        raise(ValueError("len(c2) != 3*len(a2)"))
    ## end if

    #** Confirm matching lengths of coords and atoms w/corresponding
    #  objects within the ORCA_ENGRAD instance
    if not in_g1.shape[0] == in_g2.shape[0]:
        match = False
        fail_type = "coord_dim_mismatch"
        return match, fail_type
    ## end if
    if not in_a1.shape[0] == in_a2.shape[0]:
        match = False
        fail_type = "atom_dim_mismatch"
        return match, fail_type
    ## end if

    #** Element-wise check for geometry match to within 'tol'
    fail_loc = np.less_equal(np.abs( \
                        np.subtract(in_g1,in_g2)), tol)
    if sum(fail_loc) != in_g2.shape[0]:
        # Count of matching coordinates should equal the number of
        #  coordinates. If not, complain with 'coord_mismatch' fail type.
        match = False
        fail_type = "coord_mismatch"
        return match, fail_type, fail_loc
    ## end if

    #** Element-wise check for atoms match. Quietly convert both input and
    #  instance atom arrays to atomNums to allow np.equals comparison.
    if np.issubdtype(in_a1.dtype, np.dtype('str')):
        # Presume atomic symbol data and attempt conversion
        in_a1 = np.matrix( \
                    [atomNum[e] for e in np.array(in_a1)[:,0]] \
                    ).transpose()
    ## end if
    if np.issubdtype(in_a2.dtype, np.dtype('str')):
        # Presume atomic symbol data and attempt conversion
        in_a2 = np.matrix( \
                    [atomNum[e] for e in np.array(in_a2)[:,0]] \
                    ).transpose()
    ## end if
    fail_loc = np.equal(in_a1, in_a2)

    #** Perform the test to ensure all atoms match.
    if sum(fail_loc) != in_a2.shape[0]:
        # Count of matching atoms should equal number of atoms. If not,
        #  complain with the 'atom_mismatch' fail type.
        match = False
        fail_type = "atom_mismatch"
        return match, fail_type, fail_loc

    #** If reached here, all tests passed; return success.
    return match

## end def check_geom


def template_subst(template, subs, subs_delims=['<', '>']):
    """ Perform substitution of content into tagged input template.

    No checks for valid ORCA input syntax are performed once all substitutions
    have been completed.

    Parameters
    ----------
    template : str
        Template for an ORCA input containing tags delimited by 'subs_delims',
        with tag names and substitution contents provided in 'subs'.
    subs    : iterable of two-tuples of strings
        Iterable contains tag names and corresponding content to be
        substituted into the provided template.
    subs_delims : iterable of str
        Iterable containing the 'open' and 'close' strings used to mark tags
        in the template, which are drawn from elements zero and one,
        respectively. Any elements beyond these are ignored.

    Returns
    -------
    input_text : str
        Input file generated from the parsed template, with all tag
        substitutions performed.

    Raises
    ------
    (none yet)

    """

    # Store the template into the working variable
    input_text = template

    # Iterate over subs and perform the .replace() calls
    for repset in subs:
        input_text = input_text.replace( \
                subs_delims[0] + repset[0] + subs_delims[1], repset[1])
    ##next repset

    # Return the result
    return input_text

## end def template_subst


if __name__ == '__main__':
    print("Module not executable.")





