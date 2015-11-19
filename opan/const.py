#-------------------------------------------------------------------------------
# Name:        const
# Purpose:     Module containing constants used in the various opan modules
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     13 Aug 2014
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


"""Defines objects bearing assorted constants for OpenAnharmonic.

Module-Level Members
=========================

Attributes
----------
infty : str
    Infinity symbol as Unicode string
atomNum : dict
    Atomic number lookup from element symbol

    .. note:: Keys for `atomNum` are **all uppercase** (e.g., 'AR' for argon)

atomSym : dict
    Element symbol lookup from atomic number, returned as **all uppercase**


Classes
=============

Overview
~~~~~~~~~~

Constants Classes
------------------
:class:`~opan.const.CIC` -- Application-internal code information constants

:class:`~opan.const.PHYS` -- Physical constants

:class:`~opan.const.DEF` -- Default values for parameters intended to be
user-adjustable

:class:`~opan.const.PRM` -- Internal computation parameters, intended to be
non-user-adjustable

.. :class:`~opan.const.SYMM` -- Constants relating to the point-group detection
   implementation in :mod:`opan.utils.symm`

:class:`~opan.const.UNINIT` -- Constants representing un-initialized values

:class:`~opan.const.UNITS` -- Functions returning text strings of
units descriptions

Enumeration Classes
----------------------
:class:`~opan.const.OPANEnum` -- Superclass for enumerations

    **Plain Enumerations**

    :class:`~opan.const.E_DispDirection` -- Displacement direction along
    a particular mode

    :class:`~opan.const.E_FileType` -- Various file types relevant to
    the software packages

    :class:`~opan.const.E_MassPertType`  -- Type of atomic mass
    perturbation being applied

    :class:`~opan.const.E_Software` -- Implemented computational
    software packages

    :class:`~opan.const.E_TopType` -- Molecular top classification

    **Anharmonic (VPT2) HDF5 Repository Enumerations**

    :class:`~opan.const.ERA_Data` -- Displacement-specific values

    :class:`~opan.const.ERA_Param` -- Displacement-nonspecific values

    **Units Enumerations**

    Units implemented for numerical conversions for various physical quantities.

    :class:`~opan.const.EU_RotConst` -- Rotational constants


Full Classes API
~~~~~~~~~~~~~~~~~~~

"""


# Infinity symbol
infty = u"\u221E"


class OPANEnum(object):
    """ Superclass for enumeration objects.

    Metaclassed to allow direct iteration and membership testing
    of enumeration values on the subclass type.


    .. class:: __metaclass__(type)

        Metaclass providing ability to iterate over enum values.

        With this metaclass, iterating over the class itself (rather than an
        instance) yields the valid enumeration values.

        .. method:: __iter__()

            Iterate over all defined enumeration values.

            Generator iterating over all class variables whose names match
            their contents. For a properly constructed
            :class:`~opan.error.OPANEnum` subclass, these are identical to
            the enumeration values.

            **Example:**

            >>> 'NoDisp' in E_DispDirection
            True

    """

    class __metaclass__(type):  # pragma: no cover
        def __iter__(self):
            for item in self.__dict__:
                if item == self.__dict__[item]:
                    yield item
                ## end if
            ## next item
        ## end def __iter__
    ## end class __metaclass__

## end class OPANEnum


class E_DispDirection(OPANEnum):
    """ Enumeration class for displacement directions.

    Contains enumeration parameters to indicate the displacement of the
    molecular geometry associated with a gradient, hessian or other object.

    **Enum Values**

    """
    #: Positive displacement along a particular normal mode
    Positive = 'Positive'

    #: Non-displaced geometry
    NoDisp = 'NoDisp'

    #: Negative displacement along a particular normal mode
    Negative = 'Negative'

## end class E_DispDirection


class E_MassPertType(OPANEnum):
    """ Enumeration class for atom mass perturbation types.

    Contains enumeration parameters to indicate the type of mass perturbation
    to be applied to the various atoms of a geometry, in order to: (a) break
    inertial degeneracy sufficiently to allow VPT2 computation using a lower-
    symmetry formalism; and (b), in the case of linear molecules, introduce
    an artificial 'directional preference' to the masses of the atoms
    (:attr:`ByCoord` enum value) to break the intrinsic degeneracy
    of the bending modes.

    **Enum Values**

    """

    #: Atomic masses are used without modification
    NoPerturb = 'NoPerturb'

    #: Atomic masses are perturbed atom-by-atom in an isotropic fashion
    ByAtom = 'ByAtom'

    #: Atomic masses are perturbed anisotropically, where the
    #: perturbation factor for each atom's mass varies slightly in the
    #: x-, y-, and z-directions
    ByCoord = 'ByCoord'

## end class E_MassPertType


class E_TopType(OPANEnum):
    """ Enumeration class for classifying types of molecular tops.

    Contains enumeration parameters to indicate the type of molecular top
    associated with a particular geometry.

    Inertial moments with magnitudes less than
    :attr:`opan.const.PRM.Zero_Moment_Tol` are taken as zero. Nonzero moments
    by this metric are considered to be equal if their ratio differs from
    unity by less than :attr:`opan.const.PRM.Equal_Moment_Tol`.  See
    :func:`opan.utils.inertia.principals` and the other functions defined
    in :mod:`opan.utils.inertia` for more details.

    **Enum Values**

    """

    #: Three zero principal inertial moments
    Atom = 'Atom'

    #: One zero and two equal non-zero moments
    Linear = 'Linear'

    #: Three equal, non-zero moments
    Spherical = 'Spherical'

    #: Three non-zero moments; largest two equal
    SymmProlate = 'SymmProlate'

    #: Three non-zero moments; smallest two equal
    SymmOblate = 'SymmOblate'

    #: Three unique non-zero moments
    Asymmetrical = 'Asymmetrical'

## end class E_TopType


class E_Software(OPANEnum):
    """ Enumeration class for identifying computational chemistry packages.

    This enum will be expanded if/when support for additional packages is
    implemented.

    **Enum Values**

    """

    #: The `ORCA <http://orcaforum.cec.mpg.de>`_ |external link| program package
    ORCA = 'ORCA'

## end class E_Software


class E_FileType(OPANEnum):
    """ Enumeration class for the file types generated by computational codes.

    **Enum Values**

    """

    #: XYZ atomic coordinates, assumed to follow the `Open Babel XYZ
    #: specification <http://openbabel.org/docs/2.3.0/FileFormats/
    #: XYZ_cartesian_coordinates_format.html>`_ |external link|
    xyz = 'xyz'

    #: Files containing nuclear gradient information
    grad = 'grad'

    #: Files containing nuclear Hessian information
    hess = 'hess'

    #: Files containing computational output
    output = 'output'

    #: Input files for defining computations
    inputfile = 'inputfile'

## end class E_FileType


class ERA_Data(OPANEnum):
    """ Enumeration class for datatypes in VPT2 HDF5 repository.

    Contains enumeration parameters to indicate the type of data to be
    retrieved from the on-disk HDF5 repository in the VPT2 anharmonic
    calculations of the :mod:`opan.anharm` submodule.

    **Enum Values**

    """

    #: Energy value at the given displacement
    en = 'en'

    #: Geometry at the given displacement (max precision available)
    geom = 'geom'

    #: Cartesian gradient vector
    grad = 'grad'

    #: Cartesian Hessian matrix
    hess = 'hess'

## end class ERA_Data


class ERA_Param(OPANEnum):
    """ Enumeration class for parameters in VPT2 HDF5 repository.

    Contains enumeration parameters to indicate the parameter value to be
    retrieved from the on-disk HDF5 repository in the VPT2 anharmonic
    calculations of the :mod:`opan.anharm` submodule.

    **Enum Values**

    """

    #: Length-`N` vector of all-caps atomic symbols
    atoms = 'atoms'

    #: Displacement increment in :math:`\mathrm{B}\ \mathrm{u^{1/2}}`.
    #: Note that these are *not* atomic units, which would instead be
    #: :math:`\mathrm{B}\ \mathrm{m_e^{1/2}}`.
    increment = 'increment'

    #: Cartesian center of mass of system, *with perturbations applied*,
    #: if any
    ctr_mass = 'ctr_mass'

    #: Reference values of atomic masses (unperturbed)
    ref_masses = 'ref_masses'

    #: :class:`opan.const.E_MassPertType` indicating perturbation type
    pert_mode = 'pert_mode'

    #: Length-`3*N` vector of perturbation factors (all should be ~1)
    pert_vec = 'pert_vec'

## end class ERA_Param


class EU_RotConst(OPANEnum):
    """ Units Enumeration class for rotational constants.

    Contains enumeration parameters to indicate the associated/desired units
    of interpretation/display of a rotational constant.

    String versions of these units are provided in
    :attr:`UNITS.rotConst`.

    #DOC: Add link to exposition(?) of how RotConst expression is developed,
    once written.

    **Enum Values**

    """

    #: Inverse moment of inertia, :math:`\frac{1}{\mathrm{u\ B^2}}`. Note that
    #: the mass units here are *not* atomic units, which would require
    #: :math:`\frac{1}{\mathrm{m_e\ B^2}}`.
    InvInertia = 'InvInertia'

    #: Angular frequency in atomic units, :math:`\frac{1}{\mathrm{T_a}}`
    AngFreqAtomic = 'AngFreqAtomic'

    #: Angular frequency in SI units, :math:`\frac{1}{\mathrm s}`
    #: (**NOT** :math:`\mathrm{Hz}`!)
    AngFreqSeconds = 'AngFreqSeconds'

    #: Cyclic frequency in atomic units,
    #: :math:`\frac{\mathrm{cyc}}{\mathrm{T_a}}`
    CyclicFreqAtomic = 'CyclicFreqAtomic'

    #: Cyclic frequency in :math:`\mathrm{Hz}`,
    #: :math:`\frac{\mathrm{cyc}}{\mathrm s}`
    CyclicFreqHz = 'CyclicFreqHz'

    #: Cyclic frequency in :math:`\mathrm{MHz}`,
    #: millions of :math:`\frac{\mathrm{cyc}}{\mathrm s}`
    CyclicFreqMHz = 'CyclicFreqMHz'

    #: Wavenumbers in atomic units, :math:`\frac{\mathrm{cyc}}{\mathrm{B}}`
    WaveNumAtomic = 'WaveNumAtomic'

    #: Wavenumbers in conventional units,
    #: :math:`\frac{\mathrm{cyc}}{\mathrm{cm}}`
    WaveNumCM = 'WaveNumCM'

## end class EU_RotConst


class CIC(object):
    """Container for application-internal code information constants

    These may need expansion into dictionaries keyed by
    :class:`~opan.const.E_Software` enum values, depending on the atoms
    supported by various software packages. They may also require adjustment
    to accommodate 'special' atom types such as ghost atoms and point charges.

    **Members**

    """

    #: Maximum atomic number supported
    Max_Atomic_Num = 103

    #: Minimum atomic number supported
    Min_Atomic_Num = 1

## end class CIC


class PHYS(object):

    """Container for physical constants

    **Members**

    """

    # Imports
    import numpy as _np

    #: Angstroms per Bohr radius (source: `NIST <http://physics.nist.gov/
    #: cgi-bin/cuu/Value?bohrrada0|search_for=bohr+radius>`__ |external link|)
    Ang_per_Bohr = 0.52917721067

    #: Electron mass per unified atomic mass unit (source: `NIST
    #: <http://physics.nist.gov/cgi-bin/cuu/Value?meu|
    #: search_for=electron+mass>`__ |external link|)
    me_per_amu = 1822.8885

    #: Seconds per atomic time unit (source: `NIST <http://physics.nist.gov/
    #: cgi-bin/cuu/Value?aut|search_for=atomic+time+unit>`__ |external link|)
    sec_per_Ta = 2.4188843265e-17

    #: Speed of light in atomic units, :math:`\frac{B}{T_a}`. Calculated from
    #: the `NIST <http://physics.nist.gov/cgi-bin/cuu/Value
    #: ?c|search_for=speed+of+light>`__ |external link| value for the speed of
    #: light in vacuum, :math:`2.99792458e8\ \frac{m}{s}`, using
    #: :attr:`Ang_per_Bohr` and :attr:`sec_per_Ta` as conversion factors
    light_speed = 137.036

    #: Standard Planck constant, equal to :math:`2\pi` in atomic units of
    #: :math:`\frac{\mathrm{E_h\ T_a}}{\mathrm{cyc}}`
    Planck = 2 * _np.pi

    #: Reduced Planck constant, unity by definition in the atomic units
    #: of :math:`\mathrm{E_h\ T_a}`
    Planck_bar = 1

## end class PHYS


class DEF(object):
    """Container for default parameter values (possibly user-adjustable)
    """

    from .const import E_Software as _E_SW, E_FileType as _E_FT

    #: Relative magnitude of atomic mass perturbations
    Mass_Perturbation_Magnitude = 1e-4

    #Moment_Tol = 1e-4

    #: Acceptable deviation from Kronecker delta for orthonormality testing
    Orthonorm_Tol = 1e-8

    #: Max precision of HESS geometries (currently ORCA-specific)
    HESS_Coord_Match_Tol = 1e-6

    #: Max precision of freqs in IR spectrum block (currently ORCA-specific)
    HESS_IR_Match_Tol = 1e-2

    #: Max precision of GRAD geometries (currently ORCA-specific)
    GRAD_Coord_Match_Tol = 1e-7

    #: Max tolerable deviation between XYZ geoms (currently ORCA-specific)
    XYZ_Coord_Match_Tol = 1e-12

    # Required quality of coordinate match for symmetry detection
    Symm_Match_Tol = 1e-3

    # Tolerance for deviation in searching for neighbor axes in cubic
    # symmetry groups -- value is in radians, and equals 1.5 degrees
    Symm_Axis_Match_Tol = 0.026179939

    # Rounding atomic masses to avoid precision errors in atom matching
    Symm_AtWt_Round_Digits = 4

    # Initial order of rotational symmetry to test for (conservative)
    Symm_Match_nMax = 10

    # Maximum order of atom averaging when looking for possible axes of
    # rotational symmetry
    Symm_Avg_Max = 2

    #: Dictionary of dictionaries of file extensions for geom, gradient,
    #: and hessian files from the various softwares.
    #: Access as `File_Extensions`\ [\ `E_Software`\ ][\ `E_FileType`\ ]
    File_Extensions = {
            _E_SW.ORCA :
                { _E_FT.grad : 'engrad',
                    _E_FT.hess : 'hess',
                    _E_FT.xyz : 'xyz',
                    _E_FT.output : 'out',
                    _E_FT.inputfile : 'txt'
                    }
            }

## end class DEF


class PRM(object):
    """Container for internal computation parameters (not user-adjustable)

    **Members**

    """

    #: Minimum angle deviation (degrees) required for two vectors to be
    #: considered non-parallel
    Non_Parallel_Tol = 1e-3

    #: Vector magnitude below which a vector is considered equal to the zero
    #: vector; dimensionless or units of :math:`\mathrm{B}`
    Zero_Vec_Tol = 1e-6

    #: Trap value for aberrantly large dipole  derivative values in ORCA
    #: if dipoles are not calculated in a NUMFREQ run
    Max_Sane_DipDer = 100

    #: Minimum deviation-ratio from unity below which two principal inertial
    #: moments are considered equal
    Equal_Moment_Tol = 1e-3

    #: Threshold value below which moments are  considered equal to zero;
    #: units of :math:`\mathrm{u\ B^2}`
    Zero_Moment_Tol = 1e-3

## end class PRM


class SYMM(object):
    1 # Dummy line as first line to detach the docstring
    # Be sure to also re-doc-comment the parameters in DEF and de-comment the
    # SYMM class reference in the module docstring, in addition to completing
    # this docstring, once the symmetry tools are in workable shape.
    """ Container for constants used in symmetry determination.

    Angles are in radians.

    *UNDER DEVELOPMENT*

    Td_C2_C2 = Oh_C4_C4 = 1.570796327 (90 deg)
    Td_C2_C3 = Oh_C4_C3 = 0.955316618 (54.736 deg)
    Td_C3_C3 = Oh_C3_C3 = 1.230959417 (70.529 deg)

    Icosahedral nomenclature is Ih_Cn_xCm
        n,m are orders of principal axes
        x indicates the x-th nearest Cm axis to any given Cn axis
    Ih_C5_1C5 = 1.107147967 (63.435 deg)
    Ih_C5_2C5 = 2.034442786 (116.565 deg)

    Ih_C5_1C3 = 0.651981429 (37.356 deg)
    Ih_C5_2C3 = 1.381327787 (79.144 deg)
    Ih_C5_3C3 = 1.760264867 (100.856 deg)

    Ih_C5_1C2 = 0.557442890 (31.939 deg)
    Ih_C5_2C2 = 1.024395943 (58.694 deg)
    Ih_C5_3C2 = 1.515997158 (86.860 deg)

    Ih_C3_1C3 = 0.729350896 (41.789 deg)
    Ih_C3_2C3 = 1.230262025 (70.489 deg)
    Ih_C3_3C3 = 1.910129997 (109.442 deg)

    Ih_C3_1C2 = 0.364698782 (20.896 deg)
    Ih_C3_2C2 = 0.957758721 (54.876 deg)
    Ih_C3_3C2 = 1.592672312 (91.253 deg)

    Ih_C2_1C2 = 0.632381404 (36.233 deg)
    Ih_C2_2C2 = 1.054424138 (60.414 deg)
    Ih_C2_3C2 = 1.264759664 (72.465 deg)
    Ih_C2_4C2 = 1.467100515 (84.059 deg)

    """

    Td_C2_C2 = Oh_C4_C4 = 1.570796327
    Td_C2_C3 = Oh_C4_C3 = 0.955316618
    Td_C3_C3 = Oh_C3_C3 = 1.230959417

    Ih_C5_1C5 = 1.107147967
    Ih_C5_2C5 = 2.034442786

    Ih_C5_1C3 = 0.651981429
    Ih_C5_2C3 = 1.381327787
    Ih_C5_3C3 = 1.760264867

    Ih_C5_1C2 = 0.557442890
    Ih_C5_2C2 = 1.024395943
    Ih_C5_3C2 = 1.515997158

    Ih_C3_1C3 = 0.729350896
    Ih_C3_2C3 = 1.230262025
    Ih_C3_3C3 = 1.910129997

    Ih_C3_1C2 = 0.364698782
    Ih_C3_2C2 = 0.957758721
    Ih_C3_3C2 = 1.592672312

    Ih_C2_1C2 = 0.632381404
    Ih_C2_2C2 = 1.054424138
    Ih_C2_3C2 = 1.264759664
    Ih_C2_4C2 = 1.467100515

## end class SYMM


class UNINIT(object):
    """Container for numerical values indicating an un-initialized parent
        object.

    #TODO: (?)UNINIT: Consider deprecating in favor of a custom Exception.

    **Members**
    """

    # Empty doc comments trigger inclusion in the documentation
    #:
    Unsigned_Long = -1

    #:
    Unsigned_Double = -1.0

    #:
    Signed_Long = 1234567890

    #:
    Signed_Double = 1234567890.12345

## end class UNINIT


class UNITS(object):
    """Container for dicts providing strings describing the various display
    units available for physical quantities.

    Dictionary keys are the enum values provided in the corresponding
    ``EU_xxxx`` class in this module (:mod:`opan.const`).

    ================== ===================== =====================
      Dictionary              Enum             Physical Quantity
    ================== ===================== =====================
     :attr:`rotConst`   :attr:`EU_RotConst`   Rotational constant
    ================== ===================== =====================

    """

    #TODO: (as occurs) const.UNITS: Add dicts for other units-based enums
    #TODO: (?) If good way of prettyprinting the dict becomes available,
    #   add doc comments.

    # Imports
    from .const import EU_RotConst as _EUrc

    #:
    rotConst = {
            _EUrc.InvInertia :        "1/(amu*B^2)",
            _EUrc.AngFreqAtomic :     "1/Ta",
            _EUrc.AngFreqSeconds :    "1/s",
            _EUrc.CyclicFreqAtomic :  "cyc/Ta",
            _EUrc.CyclicFreqHz :      "cyc/s",
            _EUrc.CyclicFreqMHz :     "MHz",
            _EUrc.WaveNumAtomic :     "cyc/B",
            _EUrc.WaveNumCM :         "cyc/cm"
                }

## end class UNITS


#TODO:(?) const: do atomNum/Sym need to handle ghost atoms / charges?
# Must use str.upper() when retrieving values from atomNum dictionary
atomNum = {"H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
     "F": 9, "NE": 10, "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15,
     "S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21, "TI": 22,
     "V": 23, "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29,
     "ZN": 30, "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36,
     "RB": 37, "SR": 38, "Y": 39, "ZR": 40, "NB": 41, "MO": 42, "TC": 43,
     "RU": 44, "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, "SN": 50,
     "SB": 51, "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "LA": 57,
     "CE": 58, "PR": 59, "ND": 60, "PM": 61, "SM": 62, "EU": 63, "GD": 64,
     "TB": 65, "DY": 66, "HO": 67, "ER": 68, "TM": 69, "YB": 70, "LU": 71,
     "HF": 72, "TA": 73, "W": 74, "RE": 75, "OS": 76, "IR": 77, "PT": 78,
     "AU": 79, "HG": 80, "TL": 81, "PB": 82, "BI": 83, "PO": 84, "AT": 85,
     "RN": 86, "FR": 87, "RA": 88, "AC": 89, "TH": 90, "PA": 91, "U": 92,
     "NP": 93, "PU": 94, "AM": 95, "CM": 96, "BK": 97, "CF": 98, "ES": 99,
     "FM": 100, "MD": 101, "NO": 102, "LR": 103}

atomSym = {1: 'H', 2: 'HE', 3: 'LI', 4: 'BE', 5: 'B', 6: 'C', 7: 'N',
    8: 'O', 9: 'F', 10: 'NE', 11: 'NA', 12: 'MG', 13: 'AL', 14: 'SI',
    15: 'P', 16: 'S', 17: 'CL', 18: 'AR', 19: 'K', 20: 'CA', 21: 'SC',
    22: 'TI', 23: 'V', 24: 'CR', 25: 'MN', 26: 'FE', 27: 'CO', 28: 'NI',
    29: 'CU', 30: 'ZN', 31: 'GA', 32: 'GE', 33: 'AS', 34: 'SE', 35: 'BR',
    36: 'KR', 37: 'RB', 38: 'SR', 39: 'Y', 40: 'ZR', 41: 'NB', 42: 'MO',
    43: 'TC', 44: 'RU', 45: 'RH', 46: 'PD', 47: 'AG', 48: 'CD', 49: 'IN',
    50: 'SN', 51: 'SB', 52: 'TE', 53: 'I', 54: 'XE', 55: 'CS', 56: 'BA',
    57: 'LA', 58: 'CE', 59: 'PR', 60: 'ND', 61: 'PM', 62: 'SM', 63: 'EU',
    64: 'GD', 65: 'TB', 66: 'DY', 67: 'HO', 68: 'ER', 69: 'TM', 70: 'YB',
    71: 'LU', 72: 'HF', 73: 'TA', 74: 'W', 75: 'RE', 76: 'OS', 77: 'IR',
    78: 'PT', 79: 'AU', 80: 'HG', 81: 'TL', 82: 'PB', 83: 'BI', 84: 'PO',
    85: 'AT', 86: 'RN', 87: 'FR', 88: 'RA', 89: 'AC', 90: 'TH', 91: 'PA',
    92: 'U', 93: 'NP', 94: 'PU', 95: 'AM', 96: 'CM', 97: 'BK', 98: 'CF',
    99: 'ES', 100: 'FM', 101: 'MD', 102: 'NO', 103: 'LR'}


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")
