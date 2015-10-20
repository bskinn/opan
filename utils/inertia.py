#-------------------------------------------------------------------------------
# Name:        utils.inertia
# Purpose:     Submodule containing utility functions relating to calculations
#               of inertia-related properties
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
    if geom.shape[1] != 1:
        raise(ValueError("Geometry is not a column vector"))
    ## end if
    if masses.shape[1] != 1:
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
    centered : bool
        Indicates whether the geometry has already been shifted to the
        center of mass (True --> already centered; skip pre-centering)

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



if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

