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

"""Utility functions for Anharmonic ORCA, including execution automation.

Functions implemented here are not available as of NumPy v1.8.1 and
    SciPy v0.4.9.

orthonorm_check  -- Checks orthonormality of the column vectors of a
ortho_basis      -- Constructs a 2-D basis orthogonal to a given vector
delta_fxn        -- Generalized Kronecker delta function
pack_tups        -- Pack an arbitrary combination of iterables and non-
                        iterables into a list of tuples
execute_orca     -- Execute an ORCA computation from a substitutable template
#!DOC: Complete entries for point_xxxxx functions


"""

# Imports
from orca_const import DEF, E_DispDirection


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


def check_geom(c1, a1, c2, a2, tol=DEF.XYZ_Coord_Match_Tol):
    """ Check for consistency of two geometries and atom symbol lists

    Cartesian coordinates are considered consistent with the input
        coords if each component matches to within 'tol' (default value
        specified by orca_const.DEF.Coord_Match_Tol).  If coords or
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
        from that in the ORCA_ENGRAD instance to still be considered
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
    from orca_const import atomNum
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


def execute_orca(inp_tp, work_dir, exec_cmd, subs=None, subs_delims=("<",">"), \
            sim_name="orcarun", inp_ext="txt", out_ext="out", \
            wait_to_complete=True, \
            disp_mode=0, \
            disp_dir=E_DispDirection.NoDisp, \
            disp_mag=0.0, \
            bohrs=False):
    """Executes ORCA on a dynamically constructed input file.

    Generates an ORCA input file dynamically from information passed into the
    various arguments, performs the run, and returns with exit info and
    computation results (in some fashion...).  Any required resources (GBW,
    XYZ, etc.) MUST already be present in 'work_dir'. No check for pre-existing
    files of the same base name is made; any such will be overwritten.

    ORCA MUST be called using a wrapper script; this method does not
    implement the redirection necessary to send output from a direct ORCA call
    to a file on disk.

    If 'wait_to_complete' is True, the subprocess.call() syntax will be used
    and the function will not return until ORCA execution completes. If False,
    [indicate what will be returned if not waiting.]

    #!DOC: execute_orca: The different output modes, depending on waiting or not.

    The command to call ORCA must be specified in the parameter list syntax of
    the subprocess.Popen constructor (https://docs.python.org/2/library/
    subprocess.html#popen-constructor).  The execution is flexible and general,
    to allow interface with local scripts for, e.g., submission to a job queue
    in a shared-resource environment.  The special tags "INP" and "OUT",
    enclosed with the 'subs_delims' delimiters, will be replaced with
        sim_name+'.'+inp_ext
    and
        sim_name+'.'+out_ext
    respectively in all elements of 'exec_cmd' before executing the call.
    In the special cases of inp_ext=None, the INP tag will be replaced with
    just sim_name (no extension), and similarly if out_ext=None. Similarly, the
    tag NAME will be replaced just with 'sim_name' in all cases.

    inp_ext and out_ext must be different, to avoid collisions.

    Valid ORCA input syntax of the resulting text is NOT checked before calling
    ORCA.

    No mechanism is implemented to detect hangs of ORCA. Periodic manual
    oversight is recommended.

    Parameters
    ----------
    inp_tp  : str
        Template text for the input file to be generated.
    work_dir : str
        Path to base working directory. Must already exist and contain any
        resource files (GBW, XYZ, etc.) required for the calculation.
    exec_cmd : list of str
        Sequence of strings defining the ORCA execution call in the syntax of
        the subprocess.Popen constructor. This call must be to a local script;
        stream redirection of the forked process is not supported in this
        method.
    subs    : iterable of (str, str) tuples/lists, optional
        Substitutions to be performed in template. The first element of each
        tuple corresponds to the occurrence of one or more delimited
        substitution markers in the template to be substituted for the entire
        text of the second element. Start and end delimiters are set using
        'subs_delims'.
        For example, the tuple ("ABC", "text") would convert:
            The <ABC> is working
        to:
            The text is working
        Substitutions are performed in the order they appear in the 'subs'
        tuple. Recursive substitution as the tag parsing proceeds is thus
        feasible.
    subs_delims : 2-tuple of str, optional
        Indicates starting and ending strings to be used as delimiters for the
        substitution operations into the input template; the 0th and 1st
        elements are, respectively, the starting and ending strings.  Defaults
        are "<" and ">", respectively.  Strings of any length are acceptable.
        For example, to substitute a tag of the form {|TAG|}, the tuple
        ("{|","|}") should be passed to subs_delims.  Any elements past the
        second are ignored. No checking is performed for whether the delimiters
        are "sensible" or not.
    sim_name : str, optional
        Basename to use for the input/output/working files.
        If omitted, "orcarun" will be used.
    inp_ext : str
        Extension to be used for the input file generated.
    out_ext : str
        Extension to be used for the output file generated.

    Returns
    -------
    if wait_to_complete = True:
        oo : ORCA_OUTPUT
            Collected object containing output results

    #!DOC: Return values for execute_orca, once they exist

    Raises
    ------
    ValueError : If inp_ext and out_ext are identical.

    """

    # Imports
    import os, subprocess as sp
    from orca_output import ORCA_OUTPUT
    from orca_xyz import ORCA_XYZ
    from orca_engrad import ORCA_ENGRAD
    from orca_hess import ORCA_HESS

    # Switch to dir; default exception fine for handling case of invalid dir.
    os.chdir(work_dir)

    # Check for inp_ext identical to out_ext
    if inp_ext == out_ext:
        raise(ValueError("'inp_ext' and 'out_ext' cannot be identical."))
    ##end if

    # Build the input and output file names
    if inp_ext:
        inp_fname = sim_name + '.' + inp_ext
    else:
        inp_fname = sim_name
    ##end if
    if out_ext:
        out_fname = sim_name + '.' + out_ext
    else:
        out_fname = sim_name
    ##end if

    # Perform the replacement into the exec_cmd of the input and output
    #  filenames.
    exec_cmd_subs = [ \
                s.replace(subs_delims[0] + "INP" + subs_delims[1], inp_fname) \
            for s in exec_cmd]
    exec_cmd_subs = [ \
                s.replace(subs_delims[0] + "OUT" + subs_delims[1], out_fname) \
            for s in exec_cmd_subs]
    exec_cmd_subs = [ \
                s.replace(subs_delims[0] + "NAME" + subs_delims[1], sim_name) \
            for s in exec_cmd_subs]

    # Perform the content substitutions into the template string
    input_text = template_subst(inp_tp, subs, subs_delims=subs_delims)

    # Create and write the input file
    with open(inp_fname, 'w') as input_file:
        input_file.write(input_text)
    ##end with

    # Perform the ORCA call; collect return values as appropriate
    #!TODO: execute_orca: Implement non-waiting return
    if wait_to_complete:
        # Run ORCA
        sp.call(exec_cmd_subs, cwd=os.getcwd())

        # Bind ORCA_XXXXX objects and return. Have to address possibility of
        #  things other than the output not existing
        o_out = ORCA_OUTPUT(out_fname,'file')
        try:
            o_xyz = ORCA_XYZ(path=os.path.join(work_dir, sim_name + ".xyz"), \
                                                                bohrs=bohrs)
        except IOError:
            o_xyz = None
        ## end try
        try:
            o_trj = ORCA_XYZ(path=os.path.join(work_dir, sim_name + ".trj"), \
                                                                bohrs=bohrs)
        except IOError:
            o_trj = None
        ## end try
        try:
            o_engrad = ORCA_ENGRAD(os.path.join(work_dir, sim_name + ".engrad"))
        except IOError:
            o_engrad = None
        ## end try
        try:
            o_hess = ORCA_HESS(os.path.join(work_dir, sim_name + ".hess"))
        except IOError:
            o_hess = None
        ## end try

    else:
        raise(NotImplementedError("Background execution not yet implemented."))
    ## end if

    # Return something appropriate
    #TODO: execute_orca: Must refine this, esp for the different exec modes
    return o_out, o_xyz, o_engrad, o_hess


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


class Vector(object):
    """ Container for vector operations on points or molecular geometries.

    [Assumes molecule has already been translated to center-of-mass.]

    [Molecular geometry is a vector, in order of x1, y1, z1, x2, y2, z2, ...]

    [Will need to harmonize the matrix typing; currently things are just
        passed around as np.array for the most part.]

    #!DOC: Complete Vector class docstring, including the member functions
    #!TODO: Vector: Consider casting all outputs to np.matrix()?
    """

    # Imports (those required for defaults for method parameters
    from orca_const import DEF
    import numpy as np

    # Debug var
    _SLOW = True

    # Reduced 3-D Levi-Civita matrix
    ##levi_civita = np.array([[0,1,-1],[-1,0,1],[1,-1,0]], dtype=np.float64)

    @staticmethod
    def ortho_basis(norm_vec, ref_vec=None):
        """Generates an orthonormal basis in the plane perpendicular to norm_vec

        The orthonormal basis generated spans the plane defined with norm_vec as
            its normal vector.  The handedness of on1 and on2 is such that:

                on1 x on2 == norm_vec/||norm_vec||

        norm_vec must be expressible as a one-dimensional np.array of length 3.

        Parameters
        ----------
        norm_vec : (N) vector_like
            Any numeric object expressible as an np.array of dimension one with
            length 3.  norm_vec will be converted to np.array, then raveled and
            squeezed before use.  The orthonormal basis output will span the
            plane perpendicular to norm_vec.
        ref_vec  : (N) vector_like, optional
            Same structure as norm_vec. If specified, on1 will be the normalized
            projection of ref_vec onto the plane perpendicular to norm_vec.

        Returns
        -------
        on1 : (3,1) np.matrix
            First arbitrary vector defining the orthonormal basis in the plane
            normal to vec
        on2 : (3,1) np.matrix
            Second arbitrary vector defining the orthonormal basis in the plane
            normal to vec

        Raises
        ------
        ValueError : If norm_vec or ref_vec is not expressible as a 1-D vector
            with 3 elements
        XYZError   : If ref_vec is specified and it is insufficiently non-
            parallel with respect to norm_vec
        """
        # Imports for library functions
        import numpy as np
        from numpy.random import rand
        from scipy.linalg import norm
        from scipy import arccos
        from orca_const import PRM

        # Internal parameters
        # Magnitude of the perturbation from vec in constructing a random rv
        RAND_MAG = 0.25

        # Convert to np.float64 array, ravel, and squeeze norm_vec; then test
        #  for shape and length
        nv = np.float64(np.array(norm_vec[:]).ravel().squeeze())
        if not len(nv.shape) == 1:
            raise(ValueError("norm_vec does not reduce to a 1-D array"))
        ## end if

        if not nv.shape[0] == 3:
            raise(ValueError("norm_vec length is not three"))
        ## end if

        # Normalize nv
        nv = nv / norm(nv)

        # Test for specification of ref_vec in the function call
        if ref_vec == None:
            # ref_vec not specified.
            #
            # Generate reference vector by generation of a random perturbation
            #  vector suitably non-parallel to norm_vec
            # Generate suitable randomizer, looping as needed
            rv = np.float64(1.0 - RAND_MAG + 2 * RAND_MAG * rand(3))
        ##        while 180 - np.degrees(arccos(np.dot(nv, rv) / norm(rv))) < \
        ##                                                 PRM.Non_Parallel_Tol:
            while np.degrees(arccos(abs(np.dot(nv, rv) / norm(rv)))) < \
                                                    PRM.Non_Parallel_Tol:
                rv = np.float64(1.0 - RAND_MAG + 2 * RAND_MAG * rand(3))
            ## do loop

            # Calculate perturbed vector (element-wise multiplication) and
            #  normalize
            rv = rv * nv
            rv = rv / norm(rv)

        else:
            # ref_vec specified, go ahead and use.  Start with validity check.
            rv = np.float64(np.array(ref_vec[:]).ravel().squeeze())

            if not len(rv.shape) == 1:
                raise(ValueError("ref_vec does not reduce to a 1-D array"))
            ## end if

            if not rv.shape[0] == 3:
                raise(ValueError("ref_vec length is not three"))
            ## end if

            # Normalize rv
            rv = rv / norm(rv)

            # Check for collinearity of nv and rv
    ##        if 180 - np.abs(np.degrees(arccos(np.dot(nv, rv)))) < \
    ##                                            PRM.Non_Parallel_Tol:
            # Have to put in an extra check for a dot product greater than unity
            #  due to calculation precision problems if nv == rv or nv == -rv
            if abs(np.dot(nv, rv)) > 1:
                # Essentially equal or opposite vectors, making them too nearly
                #  parallel.
                raise(XYZError("nonprl",
                        "norm_vec and ref_vec are too nearly parallel.", ""))

            if np.degrees(arccos(abs(np.dot(nv, rv)))) < \
                                                    PRM.Non_Parallel_Tol:
                # Inequal vectors, but still too nearly parallel.
                raise(XYZError("nonprl",
                        "norm_vec and ref_vec are too nearly parallel.", ""))
            ## end if

            # rv is ok to use from here

        ## end try

        # on2 is the unit vector parallel to nv x rv
        on2 = np.cross(nv, rv)
        on2 = on2 / norm(on2)

        # on1 is on2 x nv (normalization should not be necessary here)
        on1 = np.cross(on2, nv)

        # Convert to np.matrix form
        on1 = np.matrix(on1).transpose()
        on2 = np.matrix(on2).transpose()

        # Return the spanning vectors
        return on1, on2

    @staticmethod
    def orthonorm_check(a, tol=DEF.Orthonorm_Tol, report=False):
        """Checks orthonormality of the column vectors of a.

        #!DOC: [orthonorm_check: complete verbose docstring]

        Parameters
        ----------
        a : (N, M) array_like
            2D array of column vectors (all real assumed) to be checked for
            orthonormality. Does not need to be square.  N >= M strictly
            enforced, as orthonormality is only possible with <= N vectors in
            N-space.
        tol : float, optional
            Default: specified by DEF.Orthonorm_Tol of module orca_const
            Tolerance for deviation of dot products from one or zero.
        report : bool, optional
            Default: False
            Whether to record and return vectors / vector pairs failing the
            orthonormality condition.

        Returns
        -------
        o : bool
            Indicates whether column vectors of `a` are orthonormal to within
            tolerance `tol`
        n_fail : int list
            (if report == True)
            List of indices of column vectors failing normality condition. An
            empty list is returned if all vectors are normalized.
        o_fail : (int, int) list
            (if report == True)
            List of 2-tuples of indices of column vectors failing orthogonality
            condition.  An empty list is returned if all vectors are orthogonal.

        Raises
        ------
        ValueError : If an object with any non-numeric elements is provided

        ValueError : If any object is passed that is parsed as having more than
            two dimensions when converted to a matrix
        """

        # Imports
        import numpy as np

        #!TODO? orthonorm_check Must add traps to ensure a is a single array,
        #    that it is 2D, that it's all real?

        # Initialize return variables
        orth = True
        n_fail = []
        o_fail = []

        # Coerce to float64 matrix. Should handle any objects with more than
        #  two dimensions; real and all-numeric are still not yet checked, but
        #  will probably be run-time caught if too bad an object is passed.
        a_mx = np.matrix(a,dtype=np.float64)
        a_split = np.hsplit(a_mx,a_mx.shape[1])

        # Loop over vectors and check orthonormality
        for iter1 in range(a_mx.shape[1]):
            for iter2 in range(iter1,a_mx.shape[1]):
                if not abs((a_split[iter1].T * a_split[iter2])[0,0] -
                            np.float64(delta_fxn(iter1, iter2))) <= tol:
                    orth = False
                    if report:
                        if iter1 == iter2:
                            n_fail.append(iter1)
                        else:
                            o_fail.append((iter1, iter2))

        # Return results
        if report:
            return orth, n_fail, o_fail
        else:
            return orth


    @staticmethod
    def point_displ(pt1, pt2):
        """ Calculate the displacement vector between two n-D points.

        pt1 - pt2

        #!DOC: Complete point_disp docstring

        """

        #Imports
        import numpy as np

        # Make iterable
        if not np.iterable(pt1):
            pt1 = np.float64(np.array([pt1]))
        else:
            pt1 = np.float64(np.array(pt1).squeeze())
        ## end if
        if not np.iterable(pt2):
            pt2 = np.float64(np.array([pt2]))
        else:
            pt2 = np.float64(np.array(pt2).squeeze())
        ## end if

        # Calculate the displacement vector and return
        displ = np.matrix(np.subtract(pt2, pt1)).reshape(3,1)
        return displ


    @staticmethod
    def point_dist(pt1, pt2):
        """ Calculate the Euclidean distance between two n-D points.

        |pt1 - pt2|

        #!DOC: Complete point_dist docstring

        """

        # Imports
        from scipy import linalg as spla

        dist = spla.norm(Vector.point_displ(pt1, pt2))
        return dist


    @staticmethod
    def point_rotate(pt, ax, theta):
        """ Rotate a 3-D point around a 3-D axis through the origin.

        Handedness is a counter-clockwise rotation when viewing the rotation
        axis as pointing at the observer.  Thus, in a right-handed x-y-z frame,
        a 90deg rotation of (1,0,0) around the z-axis (0,0,1) yields a point at
        (0,1,0).

        #!DOC: Complete point_rotate docstring

        Raises
        ------
        ValueError : If theta is nonscalar
        ValueError : If pt or ax are not reducible to 3-D vectors
        ValueError : If norm of ax is too small
        """

        # Imports
        import numpy as np

        # Ensure pt is reducible to 3-D vector.
        pt = Vector.make_nd_vec(pt, nd=3, t=np.float64, norm=False)

        # Calculate the rotation
        rot_pt = np.dot(Vector.mtx_rot(ax, theta, reps=1), pt)

        # Should be ready to return
        return rot_pt


    @staticmethod
    def point_reflect(pt, nv):
        """ Reflect a 3-D point through a plane intersecting the origin.

        nv defines the normal vector to the plane (needs not be normalized)

        #!DOC: Complete point_reflect docstring

        Raises
        ------
        ValueError : If pt or nv are not reducible to 3-D vectors
        ValueError : If norm of nv is too small
        """

        # Imports
        import numpy as np
        from scipy import linalg as spla

        # Ensure pt is reducible to 3-D vector
        pt = Vector.make_nd_vec(pt, nd=3, t=np.float64, norm=False)

        # Transform the point and return
        refl_pt = np.dot(Vector.mtx_refl(nv, reps=1), pt)
        return refl_pt


    @staticmethod
    def geom_reflect(g, nv):
        """ Reflection symmetry operation.

        nv is normal vector to reflection plane
        g is assumed already translated to center of mass @ origin

        #!DOC: Complete Vector.geom_reflect docstring

        """

        # Imports
        import numpy as np

        # Force g to n-vector
        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)

        # Transform the geometry and return
        refl_g = np.dot(Vector.mtx_refl(nv, reps=(g.shape[0] // 3)), g) \
                    .reshape((g.shape[0],1))
        return refl_g


    @staticmethod
    def geom_rotate(g, ax, theta):
        """ Rotation symmetry operation.

        ax is rotation axis
        g is assumed already translated to center of mass @ origin

        Sense of rotation is the same as point_rotate

        #!DOC: Complete Vector.geom_rotate docstring

        """

        # Imports
        import numpy as np

        # Force g to n-vector
        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)

        # Perform rotation and return
        rot_g = np.dot(Vector.mtx_rot(ax, theta, reps=(g.shape[0] // 3)), g) \
                    .reshape((g.shape[0],1))
        return rot_g


    @staticmethod
    def symm_op(g, ax, theta, do_refl):
        """ Perform general point symmetry operation on a geometry.

        #!DOC: Complete symm_op docstring

        """

        # Imports
        import numpy as np

        # Depend on lower functions' geometry vector coercion. Just
        #  do the rotation and, if indicated, the reflection.
        gx = Vector.geom_rotate(g, ax, theta)
        if do_refl:
            gx = Vector.geom_reflect(gx, ax)
        ## end if

        # Should be good to go
        return gx


    @staticmethod
    def geom_symm_match(g, atwts, ax, theta, do_refl):
        """ [Revised match factor calculation]

        #!DOC: Complete geom_symm_match docstring

        """

        # Imports
        import numpy as np
        from scipy import linalg as spla

        # Convert g and atwts to n-D vectors
        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)
        atwts = Vector.make_nd_vec(atwts, nd=None, t=np.float64, norm=False)

        # Ensure proper dimensionality
        if not g.shape[0] == 3 * atwts.shape[0]:
            raise(ValueError("Size of 'g' is not 3*size of 'atwts'"))
        ## end if

        # Calculate transformed geometry
        gx = Vector.symm_op(g, ax, theta, do_refl)

        # Push g to a column vector
        g = g.reshape((g.shape[0],1))

        # Augment g and gx with imaginary atomic weights
        ex_wts = atwts.repeat(3,axis=0).T.reshape((atwts.shape[0]*3,1)) * 1.j
        g = np.add(g, ex_wts)
        gx = np.add(gx, ex_wts)

##        # Define calc as the outer product of the augmented vectors
##        calc = np.dot(g.reshape((g.shape[0],1)), \
##                            np.reciprocal(gx.reshape((1,gx.shape[0]))))
##
##        # Calculate the complex magnitude of each element and take log10,
##        #  then abs again
##        calc = np.abs(np.log10(np.abs(calc)))

        # Expand g and gx as column vectors of coordinates
        calc_g = g.reshape((g.shape[0] // 3, 3))
        calc_gx = gx.reshape((gx.shape[0] // 3, 3))
##
##        # Expand each into a square matrix of identical column vectors
##        calc_g = calc_g.repeat(g.shape[0], axis=1)
##        calc_gx = gx.repeat(gx.shape[0], axis=1)

        # Calc is the absolute distance between the calc-ed values,
        #  scaled by the maximum of the individual atom distances or unity.

        # Calculate the unscaled distances
        calc = [[spla.norm(np.subtract(calc_g[i,:], calc_gx[j,:])) \
                                for j in range(calc_gx.shape[0])] \
                                for i in range(calc_g.shape[0])]

        # Calculate the scale factors
        scale_g = np.array([spla.norm(calc_g[i,:]) for i in \
                    range(calc_g.shape[0])]).reshape((calc_g.shape[0],1)) \
                    .repeat(calc_g.shape[0], axis=1)
        scale_gx = np.array([spla.norm(calc_gx[j,:]) for j in \
                    range(calc_g.shape[0])]).reshape((1,calc_gx.shape[0])) \
                    .repeat(calc_gx.shape[0], axis=0)
        scale = np.maximum(np.maximum(scale_g, scale_gx),
                    np.ones_like(scale_g, dtype=np.float64))

        # Scale calc
        calc = np.divide(calc, scale)

        # Take the minimum of each row
        mins = np.min(calc, axis=1)

        # Take the maximum of the minima for the final factor
        fac = np.max(mins)

        # Using the atomic weights for checking matching can result in 'fac'
        #  being greater than unity. Return the minimum of fac and unity.
        fac = min(fac, 1.0)
        return fac


    @staticmethod
    def geom_find_rotsymm(g, atwts, ax, improp, \
            nmax=DEF.Symm_Match_nMax, \
            tol=DEF.Symm_Match_Tol):
        """ Identify highest-order symmetry for a geometry on a given axis.

        Regular and improper axes possible.

        #!DOC: Complete geom_find_rotsymm docstring

        """

        # Imports
        import numpy as np

        # Vectorize the geometry
        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)

        # Ensure a 3-D axis vector
        ax = Vector.make_nd_vec(ax, nd=3, t=np.float64, norm=True)

        # Loop downward either until a good axis is found or nval < 1
        #  Should never traverse below n == 1 for regular rotation check;
        #  could for improper, though.
        nval = nmax + 1
        nfac = 1.0
        while nfac > tol and nval > 0:
            nval = nval - 1
            try:
                nfac = Vector.geom_symm_match(g, atwts, ax, \
                                        2*np.pi/nval, improp)

            except ZeroDivisionError as zde:
                # If it's because nval == zero, ignore. Else re-raise.
                if nval > 0:
                    raise(zde)
                ## end if
            ## end try
        ## loop

        # Should be good to return
        return nval, nfac


    @staticmethod
    def geom_check_axis(g, atwts, ax, \
            nmax=DEF.Symm_Match_nMax, \
            tol=DEF.Symm_Match_Tol):
        """ [Get max proper order and reflection for an axis]

        #!DOC: Complete geom_parse_axis docstring

        """

        # Imports
        import numpy as np

        # Store the max found rotation order of the geometry.
        order = Vector.geom_find_rotsymm(g, atwts, ax, \
                                            False, nmax, tol)[0]

        # Store the presence/absence of a reflection plane.
        refl = Vector.geom_symm_match(g, atwts, ax, 0, True) < tol

        # Return the pair of values for outside handling
        return order, refl


    @staticmethod
    def geom_find_group(g, atwts, pr_ax, mom, tt, \
            nmax=DEF.Symm_Match_nMax, \
            tol=DEF.Symm_Match_Tol, \
            dig=DEF.Symm_AtWt_Round_Digits,
            avmax=DEF.Symm_Avg_Max):
        """ [Find all(?) proper rotation axes (n > 1) and reflection planes.]

        #!DOC: Complete geom_find_axes docstring INCLUDING NEW HEADER LINE

        DEPENDS on principal axes and moments being sorted such that:
            I_A <= I_B <= I_C

        Logic flow developed using:
            1) http://symmetry.otterbein.edu/common/images/flowchart.pdf
                Accessed 6 Mar 2015 (flow chart)
            2) Largent et al. J Comp Chem 22: 1637-1642 (2012).
                doi: 10.1002/jcc.22995

        Helpful examples and descriptions of point groups from:
            1) Wilson, Decius & Cross. "Molecular Vibrations." New York:
                Dover (1980), pp 82-85.
            2) "Molecular Structures of Organic Compounds -- Symmetry of
                Molecules." Website of Prof. Dr. Stefan Immel, TU Darmstadt.
                http://http://csi.chemie.tu-darmstadt.de/ak/immel/script/
                redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/
                immel/tutorials/symmetry/index7.html. Accessed 6 Mar 2015.

        Rotational symmetry numbers defined per:
            Irikura, K. K. "Thermochemistry: Appendix B: Essential Statistical
            Thermodynamics." Table II. NIST Computational Chemistry Comparison
            & Benchmark Database. Online resource: http://cccbdb.nist.gov/
            thermo.asp. Accessed 6 Mar 2015.

        """
        #!TODO: Implement principal axes threshold checking to tell if a
        #  not-strictly spherical top is far enough from spherical to ignore
        #  looking for cubic groups.  Ugh. Doesn't find the reflection planes
        #  in NH3. Going to have to explicitly deal with top type, since axes
        #  *must* be principal axes of the molecule, and off-principal axes
        #  will definitely never be symmetry elements.
        #  If asymmetric, only do pr_ax
        #  If symmetric, do the unique pr_ax and projections of atoms and
        #   midpoints normal to that axis
        #  If spherical, do everything, since every axis is inertially valid.
        #  If linear, pretty much just checking for inversion center to tell
        #   between C*v and D*h

        # Imports
        import numpy as np, itertools as itt
        from scipy import linalg as spla
        from orca_const import PRM, E_TopType as ETT
        from itertools import combinations as nCr
        from collections import namedtuple
        from orca_error import SYMMError

        # Define the Axis class
        Axis = namedtuple('Axis', 'vector order refl')

        # First, look for linear; exploit the top type, as linear should never
        #  be mis-attributed
        if tt == ETT.Linear:
            # Check for plane of symmetry; if there, D*h; if not, C*v
            #!TODO: Once symmetry element reporting structure is established,
            #  revise here to report the molecular axis as the symmetry element.
            if Vector.geom_symm_match(g, atwts, pr_ax[:,0], 0., True) < tol:
                # Has symmetry plane; D*h
                group = "D*h"
                symm_fac = 2
                return group, symm_fac
            else:
                # No symmetry plane; C*v
                group = "C*v"
                symm_fac = 1
                return group, symm_fac
            ## end if
        ## end if

        # Then, check for an atom
        if tt == ETT.Atom:
            # Simple return
            group= "Kh"
            symm_fac = 1
            return group, symm_fac
        ## end if

        # Generally, trust that the top classification is going to be more
        #  rigorous than the symmetry identification.  Thus, Spherical
        #  will almost certainly indicate a cubic group; Symmetrical, whether
        #  oblate or prolate, will indicate either a cubic group or a non-cubic
        #  with a principal rotation axis of order > 2; and Asymmetrical leaves
        #  room for any group to be found.
        # (move much of this comment to the docstring once it's working)

        # Vectorize the geometry and atwts
        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)
        atwts = Vector.make_nd_vec(atwts, nd=None, t=np.float64, norm=False)

        # Also make coordinate-split geometry
        g_coord = g.reshape((g.shape[0] // 3, 3))

        # Handle Spherical case
        if tt == ETT.Spherical:
            # Build the list of atom midpoint axes
            ax_midpts = []
            for atwt in np.unique(atwts):
                # Retrieve the sub-geometry
                g_atwt = Vector.g_subset(g, atwts, atwt, dig)

                # Only have axes to store if more than one atom
                if g_atwt.shape[0] > 3:
                    # Reshape to grouped coordinates (row vectors)
                    g_atwt = g_atwt.reshape((g_atwt.shape[0] // 3, 3))

                    # Iterate over all unique index tuples of pairs
                    for tup in nCr(range(g_atwt.shape[0]), 2):
                        # Just vector-add the appropriate atomic
                        #  coordinates; no need to normalize.
                        ax_midpts.append(np.add(*g_atwt[tup,:]))
                    ## next tup
                ## end if more than one matched atom
            ## next atwt, to index all midpoint axes in the system

            # Convert to 2-D array
            ax_midpts = np.array(ax_midpts)

            # Know for a fact that it should be a cubic group. Start looking at
            #  atom-wise vectors until an order > 1 axis is found.
            order = i = 0
            while order < 2 and i < g_coord.shape[0]:
                # Store the axis
                ax = g_coord[i,:]

                # Only check if norm is large enough
                if spla.norm(ax) > PRM.Zero_Vec_Tol:
                    order, refl = Vector.geom_check_axis(g, atwts, ax, nmax, \
                                                                        tol)
                ## end if

                # Increment
                i += 1
            ## loop

            # At this point, check to see if nothing found (could happen, e.g.
            #  in C60 buckyball) and, if not, search midpoints between like
            #  atoms, again until an order > 1 axis is found.
            #  Otherwise, store the axis information as the initial reference.
            if order >= 2:
                # Found a good axis.  Store as Axis.
                ref_Axis = Axis(vector=ax, order=order, refl=refl)
            else:
                # No good axis found along atom positions. Search midpoints.
                i = 0
                while order < 2 and i < len(ax_midpts):
                    # Store the axis
                    ax = ax_midpts[i,:]

                    # Only check if norm is large enough
                    if spla.norm(ax) > PRM.Zero_Vec_Tol:
                        order, refl = Vector.geom_check_axis(g, atwts, ax, \
                                                                    nmax, tol)
                    ## end if

                    # Increment
                    i += 1
                ## loop

                # If nothing found here, raise exception
                if order < 2:
                    raise(SYMMError(SYMMError.notfound, \
                            "Cubic point group not found in spherical top " +
                            "molecule.", "Vector.geom_find_group()"))
                ## end if

                # Store the found vector as Axis
                ref_Axis = Axis(vector=ax, order=order, refl=refl)
            ## end if

            #!RESUME: Search for other axes depending on the order of the axis found.
            return ref_Axis

            ## end if order < 2, triggering check of atom pairs

#   Leftover from originally not trusting top type
##        # Must actually search for axes &c.
##        #
##        # Initialize the container for the principal axes
##        Axes_pr = []
##        for ax in [pr_ax[:,i] for i in range(3)]:
##            order, refl = Vector.geom_check_axis(g, atwts, ax, nmax, tol)
##            if order > 1 or refl:
##                Axes_pr.append(Axis(vector=ax, order=order, refl=refl))
##            ## end if
##        ## next ax
##        return Axes_pr
##
##        # What is the max order found?
##        # If < 3, asym or sph
##        # If >=3, sym or sph; if multiple >2 then sph definitely
##
#    Not doing it this way (brute force) any more.
##        # Initialize the axes list to the principal axes (matrix of column
##        #  vectors)
##        ax_list = pr_ax
##
##        # Vectorize the geometry
##        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)
##
##        # Break into 3-vectors
##        g_vecs = np.array(np.split(g, g.shape[0] // 3))
##
##        # Add all the atom displacements to the axes list
##        ax_list = np.column_stack((ax_list, g_vecs.T))
##
##        # In each block of atom types, add axes up to 5th-order midpoints
##        for atwt in np.unique(atwts):
##            # Retrieve the sub-geometry
##            g_atwt = Vector.g_subset(g, atwts, atwt, dig)
##
##            # Reshape to grouped coordinates (row vectors)
##            g_atwt = g_atwt.reshape((g_atwt.shape[0] // 3, 3))
##
##            # If more than one atom with the given weight, start at pairs
##            #  and go up from there
##            if g_atwt.shape[0] >= 2:
##                for grp_order in range(2, 1 + min(avmax, g_atwt.shape[0])):
##                    # Retrieve all unique index tuples for the indicated order
##                    for tup in nCr(range(g_atwt.shape[0]), grp_order):
##                        # Just vector-add the appropriate atomic coordinates.
##                        #  No need to normalize or anything.
##                        ax_list = np.column_stack((ax_list, \
##                                reduce(np.add,[g_atwt[i,:] for i in tup]).T))
##                    ## next tup
##                ## next order
##            ## end if
##        ## next atwt
##
##        # Scrub any collinear axes down to uniques
##        # Filter parallel axes
##        i = 0
##        while i < ax_list.shape[1] - 1:
##            j = i + 1
##            while j < ax_list.shape[1]:
##                # For ANY collinear axes, remove until only one remains.
##                v1 = ax_list[:,i]
##                v2 = ax_list[:,j]
##                if 1 - np.abs(np.dot(v1, v2) / spla.norm(v1) / spla.norm(v2)) \
##                                                    < PRM.Non_Parallel_Tol:
##                    # Strip the duplicate vector
##                    ax_list = np.column_stack((
##                            [ax_list[:,c] for c in \
##                                    range(ax_list.shape[1]) if c <> j]
##                                                ))
##
##                    # Decrement j so that nothing is skipped
##                    j -= 1
##
##                # Increment j
##                j += 1
##            ## loop j
##
##            # Increment i
##            i += 1
##        ## loop i
##
##        # Cull any too-small axes
##        i = 0
##        while i < ax_list.shape[1]:
##            # Store vector
##            v = ax_list[:,i]
##
##            # Check magnitude
##            if spla.norm(v) < PRM.Zero_Vec_Tol:
##                # Strip if too small of magnitude
##                ax_list = np.column_stack((
##                        [ax_list[:,c] for c in \
##                                range(ax_list.shape[1]) if c <> i]
##                                            ))
##
##                # Decrement counter to maintain position in reduced array
##                i -= 1
##            ## end if
##
##            # Increment counter
##            i +=1
##        ## loop
##
##        # Search all remaining axes for rotations and reflections
##        prop_list = []
##        for v in [ax_list[:,i] for i in range(ax_list.shape[1])]:
##            order = Vector.geom_find_rotsymm(g, atwts, v, \
##                                                False, nmax, tol)[0]
##            #print("Prin: " + str(v))
##            if order > 1:
##                # Rotational axis worth reporting is found. Check reflection
##                if Vector.geom_symm_match(g, atwts, v, 0, True) < tol:
##                    # Does have a reflection
##                    prop_list.append((v,order,True))
##                else:
##                    # No reflection
##                    prop_list.append((v,order,False))
##                ## end if
##            else:
##                # No rotation, but check for reflection
##                if Vector.geom_symm_match(g, atwts, v, 0, True) < tol:
##                    # Has a reflection; do report
##                    prop_list.append((v,1,True))
##                ## end if
##            ## end if
##        ## next v
##
##        # Then test all rotations for 2x-order impropers
##
##        # Finally test for inversion center
##
##        # Then search the point group catalog and assign



        return prop_list


    @staticmethod
    def g_subset(g, atwts, atwt, \
                digits=DEF.Symm_AtWt_Round_Digits):
        """ Extract a subset of a geometry matching a desired atom.

        #!DOC: Complete g_subset docstring

        """

        # Imports
        import numpy as np

        # Ensure g and atwts are n-D vectors
        g = Vector.make_nd_vec(g, nd=None, t=np.float64, norm=False)
        atwts = Vector.make_nd_vec(atwts, nd=None, t=np.float64, norm=False)

        # Ensure dims match (should already be checked at object creation...)
        if not (len(g) == 3*len(atwts)):
            raise(ValueError("Dim mismatch [len(g) != 3*len(ats)]."))
        ## end if

        # Pull into coordinate groups
        co = np.split(g, g.shape[0] // 3)

        # Filter by the indicated atomic weight
        cf = [c for (c,a) in zip(co, atwts) if \
                        np.round(a, digits) == np.round(atwt, digits)]

        # Expand back to single vector, if possible
        if not cf == []:
            g_sub = np.concatenate(cf, axis=0)
            g_sub = g_sub.reshape((g_sub.shape[0],1))
        else:
            g_sub = []
        ## end if

        # Return the subset
        return g_sub


    @staticmethod
    def make_nd_vec(v, nd=None, t=None, norm=False):
        """ Coerce input to np.array() and validate dimensionality.

        Ensure dimensionality 'n' if passed.
        Cast to type 't' if passed
        Normalize output if norm=True

        #!DOC: Complete make_nd_vec docstring

        """

        # Imports
        import numpy as np
        from scipy import linalg as spla

        # Reduce the input to the extent possible
        out_v = np.array(v, dtype=t).squeeze()

        # Confirm vector form
        if not len(out_v.shape) == 1:
            raise(ValueError("'v' is not reducible to a vector."))
        ## end if

        # If indicated, confirm dimensionality
        if nd:
            if not out_v.shape[0] == nd:
                raise(ValueError("'v' dimension is " + str(out_v.shape[0]) + \
                        ", not " + str(nd)))
            ## end if
        ## end if

        # Normalize, if indicated
        if norm:
            out_v = out_v / spla.norm(out_v)
        ## end if

        # Return result
        return out_v


    @staticmethod
    def mtx_refl(nv, reps=1):
        """ Generate block-diagonal reflection matrix about nv.

        reps must be >=1 and indicates the number of times the reflection
        matrix should be repeated along the block diagonal.  Typically this
        will be the number of atoms in a geometry.

        #!DOC: Complete mtx_refl docstring

        """

        # Imports
        import numpy as np
        from scipy import linalg as spla
        from orca_const import PRM

        # Ensure |nv| is large enough for confident directionality
        if spla.norm(nv) < PRM.Zero_Vec_Tol:
            raise(ValueError("Norm of 'nv' is too small."))
        ## end if

        # Ensure nv is a normalized np.float64 3-vector
        nv = Vector.make_nd_vec(nv, nd=3, t=np.float64, norm=True)

        # Ensure reps is a positive scalar integer
        if not np.isscalar(reps):
            raise(ValueError("'reps' must be scalar."))
        ## end if
        if not np.issubdtype(type(reps), int):
            raise(ValueError("'reps' must be an integer."))
        ## end if
        if not reps > 0:
            raise(ValueError("'reps' must be a positive integer."))
        ## end if

        # Initialize the single-point reflection transform matrix
        base_mtx = np.zeros(shape=(3,3), dtype=np.float64)

        # Construct the single-point transform matrix
        for i in range(3):
            for j in range(i,3):
                if i==j:
                    base_mtx[i,j] = 1 - 2*nv[i]**2
                else:
                    base_mtx[i,j] = base_mtx[j,i] = -2*nv[i]*nv[j]
                ## end if
            ## next j
        ## next i

        # Construct the block-diagonal replicated reflection matrix
        refl_mtx= spla.block_diag(*[base_mtx for i in range(reps)])

        # Return the result
        return refl_mtx


    @staticmethod
    def mtx_rot(ax, theta, reps=1):
        """ Generate block-diagonal rotation matrix about ax.

        [copy handedness from somewhere]

        #!DOC: Complete mtx_rot docstring

        """

        # Imports
        import numpy as np
        from scipy import linalg as spla
        from orca_const import PRM

        # Ensure |ax| is large enough for confident directionality
        if spla.norm(ax) < PRM.Zero_Vec_Tol:
            raise(ValueError("Norm of 'ax' is too small."))
        ## end if

        # Ensure ax is a normalized np.float64 3-vector
        ax = Vector.make_nd_vec(ax, nd=3, t=np.float64, norm=True)

        # Ensure reps is a positive scalar integer
        if not np.isscalar(reps):
            raise(ValueError("'reps' must be scalar."))
        ## end if
        if not np.issubdtype(type(reps), int):
            raise(ValueError("'reps' must be an integer."))
        ## end if
        if not reps > 0:
            raise(ValueError("'reps' must be a positive integer."))
        ## end if

        # Ensure theta is scalar
        if not np.isscalar(theta):
            raise(ValueError("'theta' must be scalar."))
        ## end if

        # Assemble the modified Levi-Civita matrix
        mod_lc = np.array([ [0, -ax[2], ax[1]],
                        [ax[2], 0, -ax[0]],
                        [-ax[1], ax[0], 0] ], dtype=np.float64)

        # Compute the outer product of the axis vector
        ax_oprod = np.dot(ax.reshape((3,1)), ax.reshape((1,3)))

        # Construct the base matrix
        #  Will need to refer to external math to explain this.
        base_mtx = np.add(
                            np.add( (1.0 - np.cos(theta)) * ax_oprod,
                                                np.cos(theta) * np.eye(3)
                                  ),
                            np.sin(theta) * mod_lc
                         )

        # Construct the block-diagonal replicated reflection matrix
        rot_mtx= spla.block_diag(*[base_mtx for i in range(reps)])

        # Return the result
        return rot_mtx



if __name__ == '__main__':
    print("Module not executable.")








