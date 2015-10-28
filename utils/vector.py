#-------------------------------------------------------------------------------
# Name:        utils.vector
# Purpose:     Module containing vector/symmetry utility functions for
#               OpenAnharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     5 Oct 2015
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


""" Submodule for vector operations on

[Functions implemented here are not available as of NumPy v1.8.1 and
    SciPy v0.4.9.]

#DOC: Complete vector module docstring

"""

# Imports (those required for defaults for method parameters
from ..const import DEF as _DEF


def ortho_basis(normal, ref_vec=None):
    """Generates an orthonormal basis in the plane perpendicular to 'normal'

    The orthonormal basis generated spans the plane defined with 'normal' as
        its normal vector.  The handedness of on1 and on2 is such that:

            on1 x on2 == normal/||normal||

    'normal' must be expressible as a one-dimensional np.array of length 3.

    Parameters
    ----------
    normal   : length-3 np.float_
        The orthonormal basis output will span the plane perpendicular
        to 'normal'.
    ref_vec  : length-3 np.float_, optional
        If ref_vec specified, on1 will be the normalized projection of ref_vec
        onto the plane perpendicular to 'normal'

    Returns
    -------
    on1 : length-3 np.float_
        First vector defining the orthonormal basis in the plane
        normal to 'normal'
    on2 : length-3 np.float_
        Second vector defining the orthonormal basis in the plane
        normal to 'normal'

    Raises
    ------
    ValueError : If normal or ref_vec is not expressible as a 1-D vector
        with 3 elements
    XYZError   : If ref_vec is specified and it is insufficiently non-
        parallel with respect to 'normal'
    """

    # Imports for library functions
    import numpy as np
    from scipy import linalg as spla
    from scipy import random as sprnd
    from ..const import PRM
    from ..error import VECTORError

    # Internal parameters
    # Magnitude of the perturbation from 'normal' in constructing a random rv
    RAND_MAG = 0.25

    # Convert 'normal' to np.array, then test for shape and length
    nv = np.asarray(normal).squeeze()
    if not len(nv.shape) == 1:
        raise(ValueError("'normal' is not a vector"))
    ## end if
    if not nv.shape[0] == 3:
        raise(ValueError("Length of 'normal' is not three"))
    ## end if

    # Normalize nv
    nv = nv / spla.norm(nv)

    # Test for specification of ref_vec in the function call
    if ref_vec == None:
        # ref_vec not specified.
        #
        # Generate reference vector by generation of a random perturbation
        #  vector suitably non-parallel to norm_vec
        # Generate suitable randomizer, looping as needed
        rv = nv
        while parallel_check(nv, rv):
            rv = np.float64(1.0 - RAND_MAG + 2 * RAND_MAG * sprnd.rand(3))
        ## do loop

        # Calculate rejection of perturbed vector on the normal, then
        #  normalize
        rv = rej(rv, nv)
        rv = rv / spla.norm(rv)

    else:
        # ref_vec specified, go ahead and use.  Start with validity check.
        rv = np.asarray(ref_vec).squeeze()

        if not len(rv.shape) == 1:
            raise(ValueError("ref_vec is not a vector"))
        ## end if
        if not rv.shape[0] == 3:
            raise(ValueError("ref_vec length is not three"))
        ## end if

        # Normalize rv
        rv = rv / spla.norm(rv)

        # Check for collinearity of nv and rv; raise error if too close
        if parallel_check(nv, rv):
            # Essentially equal or opposite vectors, making them too nearly
            #  parallel.
            raise(VECTORError(VECTORError.nonprl,
                    "'normal' and ref_vec are too nearly parallel.", ""))
        ## end if

        # rv is ok to use from here

    ## end try

    # on2 is the unit vector parallel to nv x rv
    on2 = np.cross(nv, rv)
    on2 = on2 / spla.norm(on2)

    # on1 is on2 x nv (normalization should not be necessary here, but is
    #  performed just in case)
    on1 = np.cross(on2, nv)
    on1 = on1 / spla.norm(on1)

    # Return the spanning vectors
    return on1, on2

## end def ortho_basis


def orthonorm_check(a, tol=_DEF.Orthonorm_Tol, report=False):
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
        Default: specified by _DEF.Orthonorm_Tol of module .const
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
    #    that it is 2D, that it's all real? To enforce the limits stated
    #    in the docstring?

    # Initialize return variables
    orth = True
    n_fail = []
    o_fail = []

    # Coerce to float_ matrix. Should handle any objects with more than
    #  two dimensions; real and all-numeric are still not yet checked, but
    #  will probably be run-time caught if too bad an object is passed.
    a_mx = np.matrix(a,dtype=np.float_)
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

## end def orthonorm_check


def parallel_check(vec1, vec2):
    """Checks whether two vectors are parallel OR anti-parallel

    Vectors must be of the same dimension.

    Parameters
    ----------
    vec1    : length-a np.float_
        First vector to compare
    vec2    : length-a np.float_
        Second vector to compare

    Returns
    -------
    par     : bool
        True if (anti-)parallel to within 'PRM.Non_Parallel_Tol' degrees.
        False otherwise.
    """

    # Imports
    from ..const import PRM
    import numpy as np

    # Convert to squeezed np.array
    vec1 = np.array(vec1).squeeze()
    vec2 = np.array(vec2).squeeze()

    # Initialize False
    par = False

    # Shape check
    for v,n in zip([vec1, vec2], [1, 2]):
        if not len(v.shape) == 1:
            raise(ValueError("Bad shape for vector #" + str(n)))
        ## end if
    ## next v,n
    if not vec1.shape[0] == vec2.shape[0]:
        raise(ValueError("Vector length mismatch"))
    ## end if

    # Check for (anti-)parallel character and return
    angle = vec_angle(vec1, vec2)
    if min([abs(angle), abs(angle - 180.)]) < PRM.Non_Parallel_Tol:
        par = True
    ## end if

    return par

## end def parallel_check


def proj(vec, vec_onto):
    """ Vector projection

    Calculated as vec_onto * dot(vec, vec_onto) / dot(vec_onto, vec_onto)

    Parameters
    ----------
    vec      : length-N np.float_
        Vector to project
    vec_onto : length-N np.float_
        Vector onto which vec is to be projected

    Returns
    -------
    proj_vec : length-N np.float_
        Projection of vec onto vec_onto
    """

    # Imports
    import numpy as np

    # Convert to np.arrays
    vec = np.asarray(vec).squeeze()
    vec_onto = np.asarray(vec_onto).squeeze()

    # Ensure vectors
    if not len(vec.shape) == 1:
        raise(ValueError("'vec' is not a vector"))
    ## end if
    if not len(vec_onto.shape) == 1:
        raise(ValueError("'vec_onto' is not a vector"))
    ## end if
    if not vec.shape[0] == vec_onto.shape[0]:
        raise(ValueError("Shape mismatch between vectors"))
    ## end if

    # Calculate the projection and return
    proj_vec = np.float_(np.asscalar(np.dot(vec.T, vec_onto))) / \
                np.float_(np.asscalar(np.dot(vec_onto.T, vec_onto))) * vec_onto
    return proj_vec

## end def proj


def rej(vec, vec_onto):
    """ Vector rejection

    Calculated as vec - proj(vec, vec_onto).

    Parameters
    ----------
    vec      : length-N np.float_
        Vector to reject
    vec_onto : length-N np.float_
        Vector onto which vec is to be rejected

    Returns
    -------
    rej_vec : length-N np.float_
        Rejection of vec onto vec_onto
    """

    # Imports
    import numpy as np

    # Convert vec to np.array. Checking/conversion of vec_onto handled by 'proj'
    vec = np.asarray(vec).squeeze()

    # Calculate and return.
    rej_vec = vec - proj(vec, vec_onto)
    return rej_vec

## end def rej


def vec_angle(vec1, vec2):
    """ Angle between two N-dimensional vectors.

    Angle calculated as arccos(dot(vec1, vec2) / norm(vec1) / norm(vec2)).
    Should be valid for vectors in any dimension as long as they're equal.
    Any dimension mismatch raises an error.

    Parameters
    ----------
    vec1    : length-N np.array of np.float_
        First vector
    vec2    : length-N np.array of np.float_
        Second vector

    Returns
    -------
    angle   : np.float_
        Angle between the two vectors in degrees

    """

    # Imports
    import numpy as np
    from scipy import linalg as spla
    from ..const import PRM

    # Convert to np.arrays
    vec1 = np.asarray(vec1).squeeze()
    vec2 = np.asarray(vec2).squeeze()

    # Check shape and equal length
    if len(vec1.shape) != 1:
        raise(ValueError("'vec1' is not a vector"))
    ## end if
    if len(vec2.shape) != 1:
        raise(ValueError("'vec2' is not a vector"))
    ## end if
    if vec1.shape[0] != vec2.shape[0]:
        raise(ValueError("Vector lengths are not equal"))
    ## end if

    # Check magnitudes
    if spla.norm(vec1) < PRM.Zero_Vec_Tol:
        raise(ValueError("'vec1' norm is too small"))
    ## end if
    if spla.norm(vec2) < PRM.Zero_Vec_Tol:
        raise(ValueError("'vec2' norm is too small"))
    ## end if

    # Calculate the angle and return. Do in multiple steps to test for
    #  possible >1 or <-1 values from numerical precision errors.
    dotp = np.dot(vec1, vec2) / spla.norm(vec1) / spla.norm(vec2)

    if dotp > 1:
        angle = 0.
    elif dotp < -1:
        angle = 180.
    else:
        angle = np.degrees(np.arccos(dotp))
    ## end if

    return angle

## end def vec_angle


if __name__ == '__main__':
    print("Module not executable.")
