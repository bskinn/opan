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


""" Submodule for miscellaneous vector operations

Functions implemented here are (to the best of this author's knowledge)
not available in NumPy or SciPy.

**Functions**

.. autofunction:: opan.utils.vector.ortho_basis(normal[, ref_vec])

.. autofunction:: opan.utils.vector.orthonorm_check(a[, tol[, report]])

.. autofunction:: opan.utils.vector.parallel_check(vec1, vec2)

.. autofunction:: opan.utils.vector.proj(vec, vec_onto)

.. autofunction:: opan.utils.vector.rej(vec, vec_onto)

.. autofunction:: opan.utils.vector.vec_angle(vec1, vec2)

"""

# Imports (those required for defaults for method parameters
from ..const import DEF as _DEF
from .decorate import arraysqueeze as _arraysqueeze


@_arraysqueeze(0,1)
def ortho_basis(normal, ref_vec=None):
    """Generates an orthonormal basis in the plane perpendicular to `normal`

    The orthonormal basis generated spans the plane defined with `normal` as
    its normal vector.  The handedness of `on1` and `on2` in the returned
    basis is such that:

    .. math::

            \\mathsf{on1} \\times \\mathsf{on2} =
            {\\mathsf{normal} \\over \\left\\| \\mathsf{normal}\\right\\|}

    `normal` must be expressible as a one-dimensional np.array of length 3.

    Parameters
    ----------
    normal   : length-3 |npfloat|_
        The orthonormal basis output will span the plane perpendicular
        to `normal`.
    ref_vec  : length-3 |npfloat|_\\ , optional
        If specified, `on1` will be the normalized projection of `ref_vec`
        onto the plane perpendicular to `normal`. Default is |None|.

    Returns
    -------
    on1 : length-3 |npfloat|_
        First vector defining the orthonormal basis in the plane
        normal to `normal`
    on2 : length-3 |npfloat|_
        Second vector defining the orthonormal basis in the plane
        normal to `normal`

    Raises
    ------
    ~exceptions.ValueError
        If `normal` or `ref_vec` is not expressible as a 1-D vector
        with 3 elements
    ~opan.error.XYZError
        If `ref_vec` is specified and it is insufficiently non-
        parallel with respect to `normal`
    """

    # Imports for library functions
    import numpy as np
    from scipy import linalg as spla
    from scipy import random as sprnd
    from ..const import PRM
    from ..error import VectorError

    # Internal parameters
    # Magnitude of the perturbation from 'normal' in constructing a random rv
    RAND_MAG = 0.25

    # Test 'normal' for shape and length
    if not len(normal.shape) == 1:
        raise(ValueError("'normal' is not a vector"))
    ## end if
    if not normal.shape[0] == 3:
        raise(ValueError("Length of 'normal' is not three"))
    ## end if

    # Normalize to concise variable 'nv'
    nv = normal / spla.norm(normal)

    # Test for specification of ref_vec in the function call
    if ref_vec is None:
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
        if not len(ref_vec.shape) == 1:
            raise(ValueError("ref_vec is not a vector"))
        ## end if
        if not ref_vec.shape[0] == 3:
            raise(ValueError("ref_vec length is not three"))
        ## end if

        # Normalize ref_vec to 'rv'
        rv = ref_vec / spla.norm(ref_vec)

        # Check for collinearity of nv and rv; raise error if too close
        if parallel_check(nv, rv):
            # Essentially equal or opposite vectors, making them too nearly
            #  parallel.
            raise(VectorError(VectorError.NONPRL,
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


def orthonorm_check(a, tol=_DEF.ORTHONORM_TOL, report=False):
    """Checks orthonormality of the column vectors of a matrix.

    If a one-dimensional np.array is passed to `a`, it is treated as a single
    column vector, rather than a row matrix of length-one column vectors.

    The matrix `a` does not need to be square, though it must have at least
    as many rows as columns, since orthonormality is only possible in N-space
    with a set of no more than N vectors. (This condition is not directly
    checked.)

    Parameters
    ----------
    a : R x S |npfloat|_
        2-D array of column vectors to be checked for orthonormality.

    tol : |npfloat|_\\ , optional
        Tolerance for deviation of dot products from one or zero. Default
        value is :data:`opan.const.DEF.ORTHONORM_TOL`.

    report : bool, optional
        Whether to record and return vectors / vector pairs failing the
        orthonormality condition. Default is |False|.

    Returns
    -------
    o : bool
        Indicates whether column vectors of `a` are orthonormal to within
        tolerance `tol`.

    n_fail : list of |int|_ or |None|
        If `report` == |True|: A list of indices of column vectors
        failing the normality condition, or an empty list if all vectors
        are normalized.

        If `report` == |False|: |None|

    o_fail : list of 2-tuples of |int|_ or |None|
        If `report` == |True|: A list of 2-tuples of indices of
        column vectors failing the orthogonality condition, or an
        empty list if all vectors are orthogonal.

        If `report` == |False|: |None|

    """

    # Imports
    import numpy as np
    from .base import delta_fxn

    #!TODO? orthonorm_check Must add traps to ensure a is a single array,
    #    that it is 2D, that it's all real? To enforce the limits stated
    #    in the docstring?

    # Initialize return variables
    orth = True
    n_fail = []
    o_fail = []

    # Coerce to float_ matrix. Must treat 1-D vector as column vector.
    #  Should raise an exception for any objects with more than
    #  two dimensions; real and all-numeric are still not yet checked, but
    #  will probably be run-time caught if too bad an object is passed.
    if len(a.shape) == 1:
        a_mx = np.matrix(a, dtype=np.float_).T
    else:
        a_mx = np.matrix(a, dtype=np.float_)

    # Split matrix into separate vectors for convenient indexing.
    a_split = np.hsplit(a_mx, a_mx.shape[1])

    # Loop over vectors and check orthonormality.
    for iter1 in range(a_mx.shape[1]):
        for iter2 in range(iter1,a_mx.shape[1]):
            if not abs((a_split[iter1].T * a_split[iter2])[0,0] -
                        np.float_(delta_fxn(iter1, iter2))) <= tol:
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
        return orth, None, None

## end def orthonorm_check


@_arraysqueeze(0,1)
def parallel_check(vec1, vec2):
    """Checks whether two vectors are parallel OR anti-parallel.

    Vectors must be of the same dimension.

    Parameters
    ----------
    vec1    : length-R |npfloat|_
        First vector to compare
    vec2    : length-R |npfloat|_
        Second vector to compare

    Returns
    -------
    bool
        |True| if (anti-)parallel to within
        :data:`opan.const.PRM.NON_PARALLEL_TOL` degrees.  |False| otherwise.

    """

    # Imports
    from ..const import PRM
    import numpy as np

    # Initialize False
    par = False

    # Shape check
    for n,v in enumerate([vec1, vec2]):
        if not len(v.shape) == 1:
            raise(ValueError("Bad shape for vector #{0}".format(n)))
        ## end if
    ## next v,n
    if not vec1.shape[0] == vec2.shape[0]:
        raise(ValueError("Vector length mismatch"))
    ## end if

    # Check for (anti-)parallel character and return
    angle = vec_angle(vec1, vec2)
    if min([abs(angle), abs(angle - 180.)]) < PRM.NON_PARALLEL_TOL:
        par = True
    ## end if

    return par

## end def parallel_check


@_arraysqueeze(0,1)
def proj(vec, vec_onto):
    """ Vector projection.

    Calculated as:

    .. math::

         \\mathsf{vec\\_onto} * \\frac{\\mathsf{vec}\\cdot\\mathsf{vec\\_onto}}
         {\\mathsf{vec\\_onto}\\cdot\\mathsf{vec\\_onto}}

    Parameters
    ----------
    vec      : length-R |npfloat|_
        Vector to project
    vec_onto : length-R |npfloat|_
        Vector onto which `vec` is to be projected

    Returns
    -------
    length-R |npfloat|_
        Projection of `vec` onto `vec_onto`
    """

    # Imports
    import numpy as np

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


@_arraysqueeze(0,1)
def rej(vec, vec_onto):
    """ Vector rejection.

    Calculated by subtracting from `vec` the projection of `vec` onto
    `vec_onto`:

    .. math::

        \\mathsf{vec} - \\mathrm{proj}\\left(\\mathsf{vec},
        \\ \\mathsf{vec\\_onto}\\right)

    Parameters
    ----------
    vec      : length-R |npfloat|_
        Vector to reject
    vec_onto : length-R |npfloat|_
        Vector onto which `vec` is to be rejected

    Returns
    -------
    length-R |npfloat|_
        Rejection of `vec` onto `vec_onto`
    """

    # Imports
    import numpy as np

    # Calculate and return.
    rej_vec = vec - proj(vec, vec_onto)
    return rej_vec

## end def rej


@_arraysqueeze(0,1)
def vec_angle(vec1, vec2):
    """ Angle between two R-dimensional vectors.

    Angle calculated as:

    .. math::

        \\arccos\\left[
        \\frac{\\mathsf{vec1}\cdot\\mathsf{vec2}}
        {\\left\\|\\mathsf{vec1}\\right\\|
            \\left\\|\\mathsf{vec2}\\right\\|}
        \\right]

    Parameters
    ----------
    vec1    : length-R |npfloat|_
        First vector
    vec2    : length-R |npfloat|_
        Second vector

    Returns
    -------
    |npfloat|_
        Angle between the two vectors in degrees

    """

    # Imports
    import numpy as np
    from scipy import linalg as spla
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
    if spla.norm(vec1) < PRM.ZERO_VEC_TOL:
        raise(ValueError("'vec1' norm is too small"))
    ## end if
    if spla.norm(vec2) < PRM.ZERO_VEC_TOL:
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
