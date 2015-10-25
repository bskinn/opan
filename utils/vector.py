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
    from ..const import PRM

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

## end def orthonorm_check


def parallel_check(vec1, vec2):
    """Checks whether two N x 1 vectors are parallel OR anti-parallel

    Vectors MUST be column vectors.

    Parameters
    ----------
    vec1    : a x 1 np.float_
        First vector to compare
    vec2    : a x 1 np.float_
        Second vector to compare

    Returns
    -------
    par     : bool
        True if (anti-)parallel to within 'PRM.Non_Parallel_Tol' degrees.
        False otherwise.
    """

    # Imports
    from ..const import PRM
    from scipy import linalg as spla
    import scipy as sp

    # Initialize True
    par = False

    # Shape check
    for v,n in zip([vec1, vec2], range(1,3)):
        if not is_col_vec(v):
            raise(ValueError("Bad shape for vector #" + str(n)))
        ## end if
    ## next v,n
    if not vec1.shape[0] == vec2.shape[0]:
        raise(ValueError("Vector length mismatch"))
    ## end if

    # Check for (anti-)parallel character and return
    angle = sp.degrees(sp.arccos(sp.dot(vec1.T,vec2) / spla.norm(vec1) /
                                                        spla.norm(vec2)))
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
    vec      : N x 1 np.float_
        Vector to project
    vec_onto : N x 1 np.float_
        Vector onto which vec is to be projected

    Returns
    -------
    proj_vec : N x 1 np.float_
        Projection of vec onto vec_onto
    """

    # Imports
    import scipy as sp

    # Ensure column vectors
    if not is_col_vec(vec):
        raise(ValueError("'vec' is not a column vector"))
    ## end if
    if not is_col_vec(vec_onto):
        raise(ValueError("'vec_onto' is not a column vector"))
    ## end if
    if not vec.shape[0] == vec_onto.shape[0]:
        raise(ValueError("Shape mismatch between vectors"))
    ## end if

    # Calculate the projection and return
    proj_vec = sp.float_(sp.asscalar(sp.dot(vec.T, vec_onto))) / \
                sp.float_(sp.asscalar(sp.dot(vec_onto.T, vec_onto))) * vec_onto
    return proj_vec

## end def proj


def rej(vec, vec_onto):
    """ Vector rejection

    Calculated as vec - proj(vec, vec_onto).

    Parameters
    ----------
    vec      : N x 1 np.float_
        Vector to reject
    vec_onto : N x 1 np.float_
        Vector onto which vec is to be rejected

    Returns
    -------
    rej_vec : N x 1 np.float_
        Rejection of vec onto vec_onto
    """

    # Calculate and return. Size &c checking handled by 'proj'
    rej_vec = vec - proj(vec, vec_onto)
    return rej_vec

## end def proj


def is_col_vec(vec):
    """ Tests for whether an object is a column vector.

    Parameters
    ----------
    vec     : Arbitrary np.array or np.matrix
        Object to be tested for 'vectorness'

    Returns
    -------
    is_col  : bool
        True if object is a column vector, False otherwise.

    Raises
    ------
    TypeError : If 'vec' lacks a 'shape' attribute
    TypeError : If 'len(vec.shape)' raises TypeError
    TypeError : If 'vec.shape[1]' raises TypeError
    Other errors might occur if a non-Numpy type happens to have a shape
        attribute
    """

    # Initialize failure
    is_col = False

    # Check member
    try:
        s = vec.shape
    except AttributeError:
        raise(TypeError("'vec' is not a NumPy array/matrix"))
    ## end try

    # Check size
    try:
        if len(s) == 2 and s[1] == 1:
            is_col = True
        ## end if
    except TypeError:
        raise(TypeError("'vec' is not a NumPy array/matrix"))
    ## end try

    # Return
    return is_col

## end def is_col_vec


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
    import scipy as sp
    from ..const import PRM

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
    if sp.linalg.norm(vec1) < PRM.Zero_Vec_Tol:
        raise(ValueError("'vec1' norm is too small"))
    ## end if
    if sp.linalg.norm(vec2) < PRM.Zero_Vec_Tol:
        raise(ValueError("'vec2' norm is too small"))
    ## end if

    # Calculate the angle and return
    angle = sp.degrees(sp.arccos(sp.dot(vec1, vec2) /
                        sp.linalg.norm(vec1) / sp.linalg.norm(vec2)))
    return angle

## end def vec_angle


if __name__ == '__main__':
    print("Module not executable.")
