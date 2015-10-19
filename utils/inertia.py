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
    (in amu) and compute the location of the center of mass.  Center of mass
    should always be calculated on the UNPERTURBED atomic masses, and thus
    only an Nx1 vector of masses is needed.

    Parameters
    ----------
    geom    : 3N x 1 np.float_
        Coordinates of the atoms
    masses  : N x 1 np.float_
        Atomic masses of the atoms

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
    if geom.shape[0] != 3*masses.shape[0]:
        raise(ValueError("Inconsistent geometry and mass vector lengths"))
    ## end if

    # Loop over the geometry and masses with an appropriate offset into the
    #  geometry
    ctr = geom.reshape(geom.shape[0] / 3, 3).T * masses / np.sum(masses)

    # Return the vector
    return ctr

## end def ctr_mass


def ctr_geom(geom, masses):
    """ Return geometry shifted to center of mass.

    Helper function to automate / encapsulate translation of a geometry to its
    center of mass. Implemented only for MASS-UNPERTURBED geometries; thus the
    required Nx1 masses vector.

    Parameters
    ----------
    geom    : 3N x 1 np.float_
        Original coordinates of the atoms
    masses  : N x 1 np.float_
        Atomic masses of the atoms

    Returns
    -------
    ctr_geom    : 3N x 1 np.float_
        Atomic coordinates after shift to center of mass

    Raises
    -------
    ValueError  : If geom & masses shapes are inconsistent

    """

    # Imports
    import numpy as np

    # Calculate the shift vector
    shift = np.repeat(ctr_mass(geom, masses), masses.shape[0], axis=1) \
                                                .T.reshape((geom.shape[0],1))

    # Shift the geometry and return
    ctr_geom = geom - shift
    return ctr_geom

## end def ctr_geom


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

