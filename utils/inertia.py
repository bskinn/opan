#-------------------------------------------------------------------------------
# Name:        utils.inertia
# Purpose:     Submodule containing utility functions dealing with calculations
#               of inertia tensor-related properties
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     18 Oct 2015
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

"""Utilities for calculation of inertia tensor & principal axes/moments

Some of these functions may have broader applicability but are housed here
due to anticipated common usage within the context of VPT2 calculations.

Functions
---------
ctr_geom        : Shift geometry relative to center of mass
ctr_mass        : Calculate the center of mass
expand_masses   : Utility for converting an N x 1 masses vector to 3N x 1
inertia_tensor  : Calculate the inertia tensor for a system
principals      : Calculate principal axes/moments and molecular top type
"""


# Module-level imports
from .decorate import arraysqueeze as _arraysqueeze
from ..const import DEF as _DEF

# Functions

@_arraysqueeze(0,1)
def ctr_mass(geom, masses):
    """Calculate the center of mass of the indicated geometry.

    Take a geometry (coordinates in BOHRS) and atoms of the indicated masses
    (in amu) and compute the location of the center of mass.

    Parameters
    ----------
    geom    : length-3N np.float_
        Coordinates of the atoms
    masses  : length-N OR length-3N np.float_
        Atomic masses of the atoms. Length-3N option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    ctr     : length-3 np.float_
        Vector location of center of mass

    Raises
    ------
    ValueError  : If geom & masses shapes are inconsistent

    """

    # Imports
    import numpy as np
    from .base import safe_cast as scast

    # Shape check
    if len(geom.shape) != 1:
        raise(ValueError("Geometry is not a vector"))
    ## end if
    if len(masses.shape) != 1:
        raise(ValueError("Masses cannot be parsed as a vector"))
    ## end if
    if not geom.shape[0] % 3 == 0:
        raise(ValueError("Geometry is not length-3N"))
    ## end if
    if geom.shape[0] != 3*masses.shape[0] and geom.shape[0] != masses.shape[0]:
        raise(ValueError("Inconsistent geometry and masses vector lengths"))
    ## end if

    # If N masses are provided, expand to 3N; if 3N, retain.
    if geom.shape[0] == 3*masses.shape[0]:
        masses = masses.repeat(3)
    ## end if

    # Calculate the mass-weighted coordinates, reshape to group by coordinate
    #  column-wise, sum each column, then divide by the sum of masses, which
    #  must further be divided by three because there are three replicates
    #  (possibly perturbed) of the mass of each atom.
    ctr = np.multiply(geom, masses).reshape((geom.shape[0]/3,3)) \
                                .sum(axis=0).squeeze() / (masses.sum() / 3)

    # Return the vector
    return ctr

## end def ctr_mass


@_arraysqueeze(0)  # masses not used directly, so not pretreated
def ctr_geom(geom, masses):
    """ Return geometry shifted to center of mass.

    Helper function to automate / encapsulate translation of a geometry to its
    center of mass.

    Parameters
    ----------
    geom    : length-3N np.float_
        Original coordinates of the atoms
    masses  : length-N OR length-3N np.float_
        Atomic masses of the atoms. Length-3N option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    ctr_geom    : length-3N np.float_
        Atomic coordinates after shift to center of mass

    Raises
    -------
    (same as ctr_mass)

    """

    # Imports
    import numpy as np

    # Calculate the shift vector. Possible bad shape of geom or masses is
    #  addressed internally by the ctr_mass call.
    shift = np.tile(ctr_mass(geom, masses), geom.shape[0] / 3)

    # Shift the geometry and return
    ctr_geom = geom - shift
    return ctr_geom

## end def ctr_geom


@_arraysqueeze(1)  # geom reassigned to ctr_geom before use, so untreated.
def inertia_tensor(geom, masses):
    """Generate the 3x3 moment-of-inertia tensor.

    Compute the 3x3 moment-of-inertia tensor for the provided geometry
    and atomic masses.  Always recenters the geometry to the center of mass.

    Reference for moment of inertia tensor (accessed 19 Oct 2015):
        http://www.chem.ox.ac.uk/teaching/Physics%20for%20Chemists/
                                Rotation/Moment%20of%20inertia.html#Int3

    Parameters
    ----------
    geom     : length-3N np.float_
        Coordinates of the atoms
    masses   : length-N OR length-3N np.float_
        Atomic masses of the atoms. Length-3N option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    tensor   : 3 x 3 np.float_
        Moment of inertia tensor for the system

    Raises
    ------
    (same as ctr_mass)

    """

    # Imports
    import numpy as np

    # Center the geometry. Takes care of any improper shapes of geom or
    #  masses via the internal call to 'ctr_mass' within the call to 'ctr_geom'
    geom = ctr_geom(geom, masses)

    # Expand the masses if required. Shape should only ever be (N,) or (3N,),
    #  else would raise an exception within the above 'ctr_geom' call
    if geom.shape[0] == 3*masses.shape[0]:
        masses = masses.repeat(3)
    ## end if

    # Initialize the tensor matrix
    tensor = np.zeros((3,3))

    # Fill the matrix
    for i in range(3):
        for j in range(i,3):
            if i == j:
                # On-diagonal element; calculate indices to include
                ind = np.concatenate([np.array(map(lambda v: v % 3,
                                        range(i+1, i+3))) + o for o in
                                        range(0,geom.shape[0],3)])

                # Calculate the tensor element
                tensor[i,i] = np.multiply(np.square(geom[ind]),
                                                        masses[ind]).sum()
            else:
                # Off-diagonal element; calculate the indices
                ind_i = np.array(range(i,geom.shape[0]+i,3))
                ind_j = np.array(range(j,geom.shape[0]+j,3))

                # Calculate the tensor element and its symmetric partner
                tensor[i,j] = np.multiply(
                        np.sqrt(np.multiply(masses[ind_i], masses[ind_j])) ,
                        np.multiply(geom[ind_i], geom[ind_j]) ).sum() * -1
                tensor[j,i] = tensor[i,j]
            ## end if
        ## next j
    ## next i

    # Return the tensor
    return tensor

## end def inertia_tensor


@_arraysqueeze(1)  # geom reassigned to ctr_geom before use, so untreated.
def principals(geom, masses, on_tol=_DEF.Orthonorm_Tol):
    """Principal axes and moments of inertia for the indicated geometry.

    Calculated by scipy.linalg.eigh, since the moment of inertia tensor
    is symmetric (real-Hermitian) by construction.  More convenient to
    compute both the axes and moments at the same time since the eigenvectors
    must be processed to ensure repeatable results.

    The principal axes (inertia tensor eigenvectors) are processed in a
    fashion to ensure repeatable, *identical* generation, including
    orientation AND directionality.
    #DOC: Add ref to exposition in webdocs once written up.

    Parameters
    ----------
    geom     : length-3N np.float_
        Coordinates of the atoms
    masses   : length-N OR length-3N np.float_
        Atomic masses of the atoms. length-3N option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    moments : length-3 np.float_
        Principal inertial moments, sorted in increasing order
        (0 <= I_A <= I_B <= I_C)
    axes    : 3 x 3 np.float_
        Principal axes, as column vectors, sorted with the principal moments
        and processed for repeatability. The axis corresponding to moments[i]
        is retrieved as axes[:,i]
    top     : E_TopType
        Detected molecular top type

    """

    # Imports
    import numpy as np
    from scipy import linalg as spla
    from ..const import PRM, E_TopType as ETT
    from ..error import INERTIAError, VECTORError
    from .vector import rej, parallel_check as prlchk
    from .vector import orthonorm_check as orthchk

    # Center the geometry. Takes care of any improper shapes of geom or
    #  masses via the internal call to 'ctr_mass' within the call to
    #  'ctr_geom'.  Will need the centered geometry eventually anyways.
    geom = ctr_geom(geom, masses)

    # Get the inertia tensor
    tensor = inertia_tensor(geom, masses)

    # Orthogonalize and store eigenvalues/-vectors. eigh documentation says it
    #  will return ordered eigenvalues.... Store eigenvalues directly to
    #  the return variable; since eigenvectors probably need work, store them
    #  to a holding variable.
    moments, vecs = spla.eigh(tensor)

    # 'fail' init for 'top
    top = None

    # Detect top type; start with error check
    if moments[0] < -PRM.Zero_Moment_Tol:  # pragma: no cover
        # Invalid moment; raise error. Should be impossible!
        raise(INERTIAError(INERTIAError.neg_moment,
                    "Negative principal inertial moment", ""))
    elif moments[0] < PRM.Zero_Moment_Tol:
        # Zero first moment. Check whether others are too
        if all(moments < PRM.Zero_Moment_Tol):
            top = ETT.Atom
        else:
            top = ETT.Linear
        ## end if
    else:
        if abs((moments[1] / moments[0]) - 1.0) < PRM.Equal_Moment_Tol:
            # Spherical or oblate symmetrical
            if abs((moments[2] / moments[1]) - 1.0) < PRM.Equal_Moment_Tol:
                top = ETT.Spherical
            else:
                top = ETT.SymmOblate
            ## end if
        else:
            # Prolate symmetrical or Asymmetric
            if abs((moments[2] / moments[1]) - 1.0) < PRM.Equal_Moment_Tol:
                top = ETT.SymmProlate
            else:
                top = ETT.Asymmetrical
            ## end if
        ## end if
    ## end if

    # Check for nothing assigned (this should never occur!)
    if top == None:  # pragma: no cover
        raise(INERTIAError(INERTIAError.top_type,
                    "Unrecognized molecular top type",""))
    ## end if

    # Initialize the axes
    axes = np.zeros((3,3))

    # Define the axes depending on the top type
    if top == ETT.Atom:
        # Just use the coordinate axes
        axes = np.identity(3, dtype=np.float_)

    elif top == ETT.Linear:
        # Zero-moment (molecular) axis always pointed toward the first atom,
        #  or the second if the first is at center-of-mass
        if spla.norm(geom[0:3]) >= PRM.Zero_Vec_Tol:
            axes[:,0] = geom[0:3] / spla.norm(geom[0:3])
        else:
            axes[:,0] = geom[3:6] / spla.norm(geom[3:6])
        ## end if

        # Second axis is the normalized rejection of the x-axis on the first
        #  axis, unless the molecule lies along the x-axis in which case it
        #  is taken as the normalized rejection of the y-axis on the first vec.
        if prlchk(axes[:,0], np.array([1.,0.,0.])):
            # Too nearly (anti-)parallel
            axes[:,1] = rej(np.array([0.,1.,0.]), axes[:,0])
        else:
            # Sufficiently non-(anti-)parallel
            axes[:,1] = rej(np.array([1.,0.,0.]), axes[:,0])
        ## end if
        axes[:,1] /= spla.norm(axes[:,1])

        # Third axis is the first crossed with the second
        axes[:,2] = np.cross(axes[:,0], axes[:,1])

    elif top == ETT.Asymmetrical:
        # Vectors should already be orthonormal; following error should
        #  never occur
        if not orthchk(vecs, tol=on_tol):  # pragma: no cover
            raise(VECTORError(VECTORError.orthonorm,
                         "'eigh' produced non-orthonormal axes", ""))
        ## end if

        # Duplicate the vectors to the axes object
        axes = vecs.copy()

        # Orient first two axes to have positive dot products with their
        #  respective first non-zero, non-orthogonal atomic displacements.
        #  Possibly fragile to some sort of highly unusual geometry.
        axes[:,0] *= np.sign(np.dot(vecs[:,0], _fadnov(vecs[:,0], geom)))
        axes[:,1] *= np.sign(np.dot(vecs[:,1], _fadnov(vecs[:,1], geom)))

        # Orient the third axis such that a3 = a1 {cross} a2
        axes[:,2] *= np.sign(np.dot(axes[:,2], np.cross(axes[:,0], axes[:,1])))

    elif top == ETT.SymmOblate:
        # First axis is taken as the normalized rejection of the first
        #  non-orthogonal atomic displacement onto the third eigenvector.
        #  All symmetric planar molecules will fail in the call to _fadnov
        #  here, and must be treated separately.
        try:
            axes[:,0] = rej(_fadnov(vecs[:,2], geom), vecs[:,2])
        except INERTIAError as IE:
            if IE.tc == INERTIAError.bad_geom:
                # Planar molecule; just use the first nonzero atomic
                #  displacement
                planar = True
                if spla.norm(geom[0:3]) < PRM.Zero_Vec_Tol:
                    axes[:,0] = geom[3:6]
                else:
                    axes[:,0] = geom[0:3]
                ## end if
            else:
                # Some other problem; re-raise
                raise
            ## end if
        else:
            # Not planar
            planar = False

        finally:
            # Either way, normalize
            axes[:,0] /= spla.norm(axes[:,0])
        ## end try

        # If planar, take the third axis as the normalized cross product
        #  of the first and second nonzero atomic displacements. If not,
        #  take the third axis directionality as that giving a positive dot
        #  product with the first non-orthogonal atomic displacement.
        if planar:
            if spla.norm(geom[0:3]) < PRM.Zero_Vec_Tol:
                # First displacement is zero
                axes[:,2] = np.cross(geom[3:6], geom[6:9])
            elif spla.norm(geom[3:6]) < PRM.Zero_Vec_Tol:
                # Second displacement is zero
                axes[:,2] = np.cross(geom[0:3], geom[6:9])
            else:
                # First and second displacements are okay
                axes[:,2] = np.cross(geom[0:3], geom[3:6])
            ## end if

            # Regardless of which calculation, normalize the vector
            axes[:,2] /= spla.norm(axes[:,2])
        else:
            # Not planar; just reorient if necessary
            axes[:,2] = vecs[:,2] * np.sign(np.dot(vecs[:,2],
                                                    _fadnov(vecs[:,2], geom)))
        ## end if

        # Second axis is the third axis crossed with the first
        axes[:,1] = np.cross(axes[:,2], axes[:,0])

    elif top == ETT.SymmProlate:
        pass

    elif top == ETT.Spherical:
        pass
    ## end if

    # Reconfirm orthonormality. Again, the error should never occur.
    if not orthchk(axes, tol=on_tol): # pragma: no cover
        raise(VECTORError(VECTORError.orthonorm,
                    "Axis conditioning broke orthonormality",""))
    ## end if

    # Return the moments, axes, and top type
    return moments, axes, top

##end def principals


@_arraysqueeze(0,1)
def _fadnov(vec, geom):
    """First non-zero Atomic Displacement Non-Orthogonal to Vec

    Utility function to identify the first atomic displacement in a geometry
    that is (a) not the zero vector; and (b) not normal to the reference vector.

    Parameters
    ----------
    vec     : length-3 np.float_
        Reference vector. Does not need to be normalized
    geom    : length-3N np.float_
        *CENTERED* molecular geometry

    Returns
    -------
    out_vec : length-3 np.float_
        Normalized non-zero atomic displacement not orthogonal to vec

    """

    # Imports
    import numpy as np
    from scipy import linalg as spla
    from ..const import PRM
    from ..error import INERTIAError
    from .vector import orthonorm_check as onchk

    # Geom and vec must both be the right shape
    if not (len(geom.shape) == 1 and geom.shape[0] % 3 == 0):
        raise(ValueError("Geometry is not length 3N"))
    ## end if
    if not vec.shape == (3,):
        raise(ValueError("Reference vector is not length 3"))
    ## end if

    # vec must not be the zero vector
    if spla.norm(vec) < PRM.Zero_Vec_Tol:
        raise(ValueError("Reference vector norm is too small"))
    ## end if

    # Normalize the ref vec
    vec = vec / spla.norm(vec)

    # Iterate over reshaped geometry
    for disp in geom.reshape((geom.shape[0]/3, 3)):
        # See if the displacement is nonzero
        if spla.norm(disp) >= PRM.Zero_Vec_Tol:
            # See if it's not orthonormal
            if not onchk(np.column_stack((disp / spla.norm(disp),
                                                vec / spla.norm(vec)))):
                # This is the displacement you are looking for
                out_vec = disp / spla.norm(disp)
                return out_vec
            ## end if
        ## end if
    ## next disp
    else:
        # Nothing fit the bill - must be atom, linear, or planar
        raise(INERTIAError(INERTIAError.bad_geom,
                    "No suitable atomic displacement found", ""))
    ## end for disp

## end def _fadnpv



if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

