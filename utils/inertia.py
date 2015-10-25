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


# Functions

def ctr_mass(geom, masses):
    """Calculate the center of mass of the indicated geometry.

    Take a geometry (coordinates in BOHRS) and atoms of the indicated masses
    (in amu) and compute the location of the center of mass.

    Parameters
    ----------
    geom    : 3N x 1 np.float_
        Coordinates of the atoms
    masses  : N x 1 OR 3N x 1 np.float_
        Atomic masses of the atoms. 3N x 1 option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    ctr     : 3 x 1 np.matrix of np.float_
        Vector location of center of mass

    Raises
    ------
    ValueError  : If geom & masses shapes are inconsistent

    """

    # Imports
    import numpy as np

    # Shape check
    if len(geom.shape) != 2 or geom.shape[1] != 1:
        raise(ValueError("Geometry is not a column vector"))
    ## end if
    if len(masses.shape) != 2 or masses.shape[1] != 1:
        raise(ValueError("Masses are not a column vector"))
    ## end if
    if not geom.shape[0] % 3 == 0:
        raise(ValueError("Geometry is not 3N x 1"))
    ## end if
    if geom.shape[0] != 3*masses.shape[0] and geom.shape[0] != masses.shape[0]:
        raise(ValueError("Inconsistent geometry and mass vector lengths"))
    ## end if

    # If N x 1 masses are provided, expand; if 3N x 1, retain.
    if geom.shape[0] == 3*masses.shape[0]:
        masses = expand_masses(masses)
    ## end if

    # Calculate the mass-weighted coordinates, reshape and transpose to group
    #  by coordinate row-wise, sum each row, then divide by the sum of masses,
    #  which must be divided by three because there are three replicates
    #  (possibly perturbed) of the mass of each atom.
    ctr = np.multiply(geom, masses).reshape((geom.shape[0]/3,3)).T \
                .sum(axis=1) / (masses.sum() / 3)

    # Return the vector
    return ctr

## end def ctr_mass


def ctr_geom(geom, masses):
    """ Return geometry shifted to center of mass.

    Helper function to automate / encapsulate translation of a geometry to its
    center of mass.

    Parameters
    ----------
    geom    : 3N x 1 np.float_
        Original coordinates of the atoms
    masses  : N x 1 OR 3N x 1 np.float_
        Atomic masses of the atoms. 3N x 1 option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    ctr_geom    : 3N x 1 np.float_
        Atomic coordinates after shift to center of mass

    Raises
    -------
    (same as ctr_mass)

    """

    # Imports
    import numpy as np

    # Calculate the shift vector. Possible bad shape of geom or masses is
    #  addressed internally by the ctr_mass call.
    shift = np.repeat(ctr_mass(geom, masses), geom.shape[0] / 3, axis=1) \
                                                .T.reshape((geom.shape[0],1))

    # Shift the geometry and return
    ctr_geom = geom - shift
    return ctr_geom

## end def ctr_geom


def inertia_tensor(geom, masses):
    """Generate the 3x3 moment-of-inertia tensor.

    Compute the 3x3 moment-of-inertia tensor for the provided geometry
    and atomic masses.  Always recenters the geometry to the center of mass.

    Reference for moment of inertia tensor (accessed 19 Oct 2015):
        http://www.chem.ox.ac.uk/teaching/Physics%20for%20Chemists/
                                Rotation/Moment%20of%20inertia.html#Int3

    Parameters
    ----------
    geom     : 3N x 1 np.float_
        Coordinates of the atoms
    masses   : N x 1 OR 3N x 1 np.float_
        Atomic masses of the atoms. 3N x 1 option is to allow calculation of
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

    # Expand the masses if required. Shape should only ever be N x 1 or 3N x 1
    if geom.shape[0] == 3*masses.shape[0]:
        masses = expand_masses(masses)
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
                tensor[i,i] = np.multiply(np.square(geom[ind,0]),
                                                        masses[ind,0]).sum()
            else:
                # Off-diagonal element; calculate the indices
                ind_i = np.array(range(i,geom.shape[0]+i,3))
                ind_j = np.array(range(j,geom.shape[0]+j,3))

                # Calculate the tensor element and its symmetric partner
                tensor[i,j] = np.multiply(
                        np.sqrt(np.multiply(masses[ind_i,0], masses[ind_j,0])) ,
                        np.multiply(geom[ind_i,0], geom[ind_j,0]) ).sum() * -1
                tensor[j,i] = tensor[i,j]
            ## end if
        ## next j
    ## next i

    # Return the tensor
    return tensor

## end def inertia_mtx


def principals(geom, masses):
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
    geom     : 3N x 1 np.float_
        Coordinates of the atoms
    masses   : N x 1 OR 3N x 1 np.float_
        Atomic masses of the atoms. 3N x 1 option is to allow calculation of
        a per-coordinate perturbed value.

    Returns
    -------
    moments : 1 x 3 np.float_
        Principal inertial moments, sorted in increasing order
        (0 <= I_A <= I_B <= I_C)
    axes    : 3 x 3 np.float)
        Principal axes, sorted with the principal moments and processed
        for repeatability
    top     : E_TopType
        Detected molecular top type

    """

    # Imports
    import numpy as np
    import scipy as sp
    from scipy import linalg as spla
    from ..const import PRM, E_TopType as ETT
    from ..error import INERTIAError
    from .vector import rej, parallel_check as prlchk

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
    if moments[0] < -PRM.Zero_Moment_Tol:
        # Invalid moment; raise error
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
        if abs(moments[0] - moments[1]) < PRM.Equal_Moment_Tol:
            # Spherical or oblate symmetrical
            if abs(moments[1] - moments[2]) < PRM.Equal_Moment_Tol:
                top = ETT.Spherical
            else:
                top = ETT.SymmOblate
            ## end if
        else:
            # Prolate symmetrical or Asymmetric
            if abs(moments[1] - moments[2]) < PRM.Equal_Moment_Tol:
                top = ETT.SymmProlate
            else:
                top = ETT.Asymmetrical
            ## end if
        ## end if
    ## end if

    # Check for nothing assigned
    if top == None:
        raise(INERTIAError(INERTIAError.top_type,
                    "Unrecognized molecular top type",""))
    ## end if

    # Initialize the axes
    axes = sp.zeros((3,3))

    # Define the axes depending on the top type
    if top == ETT.Atom:
        # Just use the coordinate axes
        axes = sp.identity(3, dtype=np.float_)

    elif top == ETT.Linear:
        # Zero-moment (molecular) axis always pointed toward the first atom,
        #  or the second if the first is at center-of-mass
        if spla.norm(geom[0:3,0]) >= PRM.Zero_Vec_Tol:
            axes[:,(0,)] = geom[0:3,0] / spla.norm(geom[0:3,0])
        else:
            axes[:,(0,)] = geom[3:6,0] / spla.norm(geom[3:6,0])
        ## end if

        # Second axis is the normalized rejection of the x-axis on the first
        #  axis, unless the molecule lies along the x-axis in which case it
        #  is taken as the normalized rejection of the y-axis on the first vec.
        #  The full rejection is calculated in the latter case to afford a 
        #  slight increase in accuracy, since the axis may deviate slightly
        #  from the x-axis
        if prlchk(axes[:,(0,)], np.array([[1.,0.,0.]]).T):
            # Too nearly (anti-)parallel
            axes[:,(1,)] = rej(np.array([[0.,1.,0.]]).T, axes[:,(0,)])
        else:
            # Sufficiently non-(anti-)parallel
            axes[:,(1,)] = rej(np.array([[1.,0.,0.]]).T, axes[:,(0,)])
        ## end if
        axes[:,1] /= spla.norm(axes[:,1])

        ## Third axis is the first crossed with the second
        axes[:,2] = sp.cross(axes[:,0], axes[:,1])

    elif top == ETT.Asymmetrical:
        pass

    elif top == ETT.SymmOblate:
        pass

    elif top == ETT.SymmProlate:
        pass

    elif top == ETT.Spherical:
        pass
    ## end if

    # Matrixify moments
    moments = np.asmatrix(moments)

    # Return the moments, axes, and top type
    return moments, axes, top

##end def principals


def expand_masses(masses):
    """ Replicate an N x 1 vector of masses to a 3N x 1 vector.

    Helper function for expanding a non-per-coordinate-perturbed masses vector
    to a full 3N x 1 dimension.

    Parameters
    ----------
    masses   : N x 1 np.float_
        Vector of masses to expand

    Returns
    -------
    expanded : 3N x 1 np.float_
        Expanded vector of masses

    """

    # Calculate the expanded vector (implicitly assumes numpy array) and
    #  return
    expanded = masses.repeat(3, axis=1).reshape((3*masses.shape[0],1))
    return expanded

## end def expand_masses


def _fadnpv(vec, geom):
    """First non-zero Atomic Displacement that is Non-Parallel with Vec

    Utility function to identify the first atomic displacement in a geometry
    that is both (a) not the zero vector and (b) non-(anti-)parallel with a
    reference vector.

    Parameters
    ----------
    vec     : 3 x 1 np.float_
        Reference vector. Does not need to be normalized
    geom    : 3N x 1 np.float_
        *CENTERED* molecular geometry

    Returns
    -------
    out_vec : 3 x 1 np.float_
        Normalized non-zero, non-atomic displacement not (anti-)parallel to vec

    """

    # Imports
    import numpy as np
    import scipy as sp
    from scipy import linalg as spla
    from ..const import PRM
    from ..error import INERTIAError
    from .vector import parallel_check as parchk

    # Geom and vec must both be the right shape
    if not (len(geom.shape) == 2 and geom.shape[0] % 3 == 0 and
                                                        geom.shape[1] == 1):
        raise(ValueError("Geometry is not 3N x 1"))
    ## end if
    if not vec.shape == (3,1):
        raise(ValueError("Reference vector is not 3 x 1"))
    ## end if

    # vec must not be the zero vector
    if spla.norm(vec) < PRM.Zero_Vec_Tol:
        raise(ValueError("Reference vector norm is too small"))
    ## end if

    # Array-ify and normalize the ref vec
    vec = np.asarray(vec).squeeze() / spla.norm(vec)

    # Iterate over reshaped geometry
    for disp in np.asarray(geom).reshape((geom.shape[0]/3, 3)):
        # See if the displacement is nonzero
        if spla.norm(disp) >= PRM.Zero_Vec_Tol:
            # See if it's nonparallel to the ref vec
            if not parchk(disp.reshape((3,1)), vec):
                # This is the displacement you are looking for
                out_vec = np.matrix(disp / spla.norm(disp)).reshape((3,1))
                break
            ## end if
        ## end if
    ## next disp
    else:
        # Nothing fit the bill - must be a linear molecule?
        raise(INERTIAError(INERTIAError.linear_mol,
                    "Linear molecule, no non-parallel displacement", ""))
    ## end for disp

    # Return the resulting vector
    return out_vec

## end def _fadnpv



if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

