#-------------------------------------------------------------------------------
# Name:        utils.symm
# Purpose:     Module containing molecular symmetry utility functions for
#               OpenAnharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     19 Oct 2015
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


""" Utilities submodule for molecular symmetry operations and detection.

[Assumes molecule has already been translated to center-of-mass.]

[Molecular geometry is a vector, in order of x1, y1, z1, x2, y2, z2, ...]

[Will need to harmonize the matrix typing; currently things are just
    passed around as np.array for the most part.]

#!DOC: Complete symm module docstring, including the member functions
#!TODO: Vector: Consider casting all outputs to np.matrix()?
"""

# Imports
from ..const import DEF as _DEF


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

## end def point_displ


def point_dist(pt1, pt2):
    """ Calculate the Euclidean distance between two n-D points.

    |pt1 - pt2|

    #!DOC: Complete point_dist docstring

    """

    # Imports
    from scipy import linalg as spla

    dist = spla.norm(point_displ(pt1, pt2))
    return dist

## end def point_dist


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
    pt = make_nd_vec(pt, nd=3, t=np.float64, norm=False)

    # Calculate the rotation
    rot_pt = np.dot(mtx_rot(ax, theta, reps=1), pt)

    # Should be ready to return
    return rot_pt

## end def point_rotate


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
    pt = make_nd_vec(pt, nd=3, t=np.float64, norm=False)

    # Transform the point and return
    refl_pt = np.dot(mtx_refl(nv, reps=1), pt)
    return refl_pt

## end def point_reflect


def geom_reflect(g, nv):
    """ Reflection symmetry operation.

    nv is normal vector to reflection plane
    g is assumed already translated to center of mass @ origin

    #!DOC: Complete geom_reflect docstring

    """

    # Imports
    import numpy as np

    # Force g to n-vector
    g = make_nd_vec(g, nd=None, t=np.float64, norm=False)

    # Transform the geometry and return
    refl_g = np.dot(mtx_refl(nv, reps=(g.shape[0] // 3)), g) \
                .reshape((g.shape[0],1))
    return refl_g

## end def geom_reflect


def geom_rotate(g, ax, theta):
    """ Rotation symmetry operation.

    ax is rotation axis
    g is assumed already translated to center of mass @ origin

    Sense of rotation is the same as point_rotate

    #!DOC: Complete geom_rotate docstring

    """

    # Imports
    import numpy as np

    # Force g to n-vector
    g = make_nd_vec(g, nd=None, t=np.float64, norm=False)

    # Perform rotation and return
    rot_g = np.dot(mtx_rot(ax, theta, reps=(g.shape[0] // 3)), g) \
                .reshape((g.shape[0],1))
    return rot_g

## end def geom_rotate


def symm_op(g, ax, theta, do_refl):
    """ Perform general point symmetry operation on a geometry.

    #!DOC: Complete symm_op docstring

    """

    # Imports
    import numpy as np

    # Depend on lower functions' geometry vector coercion. Just
    #  do the rotation and, if indicated, the reflection.
    gx = geom_rotate(g, ax, theta)
    if do_refl:
        gx = geom_reflect(gx, ax)
    ## end if

    # Should be good to go
    return gx

## end def symm_op


def geom_symm_match(g, atwts, ax, theta, do_refl):
    """ [Revised match factor calculation]

    #!DOC: Complete geom_symm_match docstring

    """

    # Imports
    import numpy as np
    from scipy import linalg as spla

    # Convert g and atwts to n-D vectors
    g = make_nd_vec(g, nd=None, t=np.float64, norm=False)
    atwts = make_nd_vec(atwts, nd=None, t=np.float64, norm=False)

    # Ensure proper dimensionality
    if not g.shape[0] == 3 * atwts.shape[0]:
        raise(ValueError("Size of 'g' is not 3*size of 'atwts'"))
    ## end if

    # Calculate transformed geometry
    gx = symm_op(g, ax, theta, do_refl)

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

## end def geom_symm_match


def geom_find_rotsymm(g, atwts, ax, improp, \
        nmax=_DEF.Symm_Match_nMax, \
        tol=_DEF.Symm_Match_Tol):
    """ Identify highest-order symmetry for a geometry on a given axis.

    Regular and improper axes possible.

    #!DOC: Complete geom_find_rotsymm docstring

    """

    # Imports
    import numpy as np

    # Vectorize the geometry
    g = make_nd_vec(g, nd=None, t=np.float64, norm=False)

    # Ensure a 3-D axis vector
    ax = make_nd_vec(ax, nd=3, t=np.float64, norm=True)

    # Loop downward either until a good axis is found or nval < 1
    #  Should never traverse below n == 1 for regular rotation check;
    #  could for improper, though.
    nval = nmax + 1
    nfac = 1.0
    while nfac > tol and nval > 0:
        nval = nval - 1
        try:
            nfac = geom_symm_match(g, atwts, ax, \
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

## end def geom_find_rotsymm


def geom_check_axis(g, atwts, ax, \
        nmax=_DEF.Symm_Match_nMax, \
        tol=_DEF.Symm_Match_Tol):
    """ [Get max proper order and reflection for an axis]

    #!DOC: Complete geom_parse_axis docstring

    """

    # Imports
    import numpy as np

    # Store the max found rotation order of the geometry.
    order = geom_find_rotsymm(g, atwts, ax, \
                                        False, nmax, tol)[0]

    # Store the presence/absence of a reflection plane.
    refl = geom_symm_match(g, atwts, ax, 0, True) < tol

    # Return the pair of values for outside handling
    return order, refl

## end def geom_check_axis


def geom_find_group(g, atwts, pr_ax, mom, tt, \
        nmax=_DEF.Symm_Match_nMax, \
        tol=_DEF.Symm_Match_Tol, \
        dig=_DEF.Symm_AtWt_Round_Digits,
        avmax=_DEF.Symm_Avg_Max):
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
    from ..const import PRM, E_TopType as ETT
    from itertools import combinations as nCr
    from collections import namedtuple
    from ..error import SYMMError

    # Define the Axis class
    Axis = namedtuple('Axis', 'vector order refl')

    # First, look for linear; exploit the top type, as linear should never
    #  be mis-attributed
    if tt == ETT.Linear:
        # Check for plane of symmetry; if there, D*h; if not, C*v
        #!TODO: Once symmetry element reporting structure is established,
        #  revise here to report the molecular axis as the symmetry element.
        if geom_symm_match(g, atwts, pr_ax[:,0], 0., True) < tol:
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
    g = make_nd_vec(g, nd=None, t=np.float64, norm=False)
    atwts = make_nd_vec(atwts, nd=None, t=np.float64, norm=False)

    # Also make coordinate-split geometry
    g_coord = g.reshape((g.shape[0] // 3, 3))

    # Handle Spherical case
    if tt == ETT.Spherical:
        # Build the list of atom midpoint axes
        ax_midpts = []
        for atwt in np.unique(atwts):
            # Retrieve the sub-geometry
            g_atwt = g_subset(g, atwts, atwt, dig)

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
                order, refl = geom_check_axis(g, atwts, ax, nmax, \
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
                    order, refl = geom_check_axis(g, atwts, ax, \
                                                                nmax, tol)
                ## end if

                # Increment
                i += 1
            ## loop

            # If nothing found here, raise exception
            if order < 2:
                raise(SYMMError(SYMMError.notfound, \
                        "Cubic point group not found in spherical top " +
                        "molecule.", "geom_find_group()"))
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
##            order, refl = geom_check_axis(g, atwts, ax, nmax, tol)
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
##        g = make_nd_vec(g, nd=None, t=np.float64, norm=False)
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
##            g_atwt = g_subset(g, atwts, atwt, dig)
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
##            order = geom_find_rotsymm(g, atwts, v, \
##                                                False, nmax, tol)[0]
##            #print("Prin: " + str(v))
##            if order > 1:
##                # Rotational axis worth reporting is found. Check reflection
##                if geom_symm_match(g, atwts, v, 0, True) < tol:
##                    # Does have a reflection
##                    prop_list.append((v,order,True))
##                else:
##                    # No reflection
##                    prop_list.append((v,order,False))
##                ## end if
##            else:
##                # No rotation, but check for reflection
##                if geom_symm_match(g, atwts, v, 0, True) < tol:
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

## end def geom_find_group


def g_subset(g, atwts, atwt, \
            digits=_DEF.Symm_AtWt_Round_Digits):
    """ Extract a subset of a geometry matching a desired atom.

    #!DOC: Complete g_subset docstring

    """

    # Imports
    import numpy as np

    # Ensure g and atwts are n-D vectors
    g = make_nd_vec(g, nd=None, t=np.float64, norm=False)
    atwts = make_nd_vec(atwts, nd=None, t=np.float64, norm=False)

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

## end def g_subset


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

## end def make_nd_vec


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
    from ..const import PRM

    # Ensure |nv| is large enough for confident directionality
    if spla.norm(nv) < PRM.Zero_Vec_Tol:
        raise(ValueError("Norm of 'nv' is too small."))
    ## end if

    # Ensure nv is a normalized np.float64 3-vector
    nv = make_nd_vec(nv, nd=3, t=np.float64, norm=True)

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

## end def mtx_refl


def mtx_rot(ax, theta, reps=1):
    """ Generate block-diagonal rotation matrix about ax.

    [copy handedness from somewhere]

    #!DOC: Complete mtx_rot docstring

    """

    # Imports
    import numpy as np
    from scipy import linalg as spla
    from ..const import PRM

    # Ensure |ax| is large enough for confident directionality
    if spla.norm(ax) < PRM.Zero_Vec_Tol:
        raise(ValueError("Norm of 'ax' is too small."))
    ## end if

    # Ensure ax is a normalized np.float64 3-vector
    ax = make_nd_vec(ax, nd=3, t=np.float64, norm=True)

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

## end def mtx_rot



if __name__ == '__main__':
    print("Module not executable.")
