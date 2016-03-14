#-------------------------------------------------------------------------------
# Name:        const
# Purpose:     Module containing constants used in the various opan modules
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     13 Aug 2014
# Copyright:   (c) Brian Skinn 2016
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
infty
    |str| -- Infinity symbol from Unicode

atom_num
    |dict| -- Atomic number lookup from element symbol

    .. note:: Keys for `atom_num` are **all uppercase** (e.g., 'AR' for argon)

atom_sym
    |dict| -- Element symbol lookup from atomic number,
    returned as **all uppercase**


Classes
=============

Overview
~~~~~~~~~~

Constants Classes
------------------
:class:`~opan.const.CIC` -- Application-internal code information constants

:class:`~opan.const.DEF` -- Default values for parameters intended to be
user-adjustable

:class:`~opan.const.PHYS` -- Physical constants

:class:`~opan.const.PRM` -- Internal computation parameters, intended to be
non-user-adjustable

.. :class:`~opan.const.SYMM` -- Constants relating to the point-group detection
   implementation in :mod:`opan.utils.symm`

:class:`~opan.const.UNINIT` -- Constants representing un-initialized values

:class:`~opan.const.UNITS` -- Functions returning text strings of
units descriptions

Enumeration Classes
----------------------
:class:`~opan.const.OpanEnum` -- Superclass for enumerations

    **Plain Enumerations**

    :class:`~opan.const.EnumDispDirection` -- Displacement direction along
    a particular mode

    :class:`~opan.const.EnumFileType` -- Various file types relevant to
    the software packages

    :class:`~opan.const.EnumMassPertType`  -- Type of atomic mass
    perturbation being applied

    :class:`~opan.const.EnumSoftware` -- Implemented computational
    software packages

    :class:`~opan.const.EnumTopType` -- Molecular top classification

    **Anharmonic (VPT2) HDF5 Repository Enumerations**

    :class:`~opan.const.EnumAnharmRepoData` -- Displacement-specific values

    :class:`~opan.const.EnumAnharmRepoParam` -- Displacement-nonspecific values

    **Units Enumerations**

    Units implemented for numerical conversions for various physical quantities.

    :class:`~opan.const.EnumUnitsRotConst` -- Rotational constants


API
~~~~

"""


# Infinity symbol
infty = "\u221E"


# ======  Enums  ====== #

class OpanEnum(object):
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
            :class:`~opan.error.OpanEnum` subclass, these are identical to
            the enumeration values.

            **Example:**

            >>> [val for val in opan.const.EnumDispDirection]
            ['POSITIVE', 'NO_DISP', 'NEGATIVE']

        .. method:: __contains__(value)

            Returns |True| if `value` is a valid value for the
            enumeration type, else |False|.

            **Example:**

            >>> 'NO_DISP' in EnumDispDirection
            True

    """

    class __metaclass__(type):
        def __iter__(self):
            for item in self.__dict__:
                if item == self.__dict__[item]:
                    yield item
                ## end if
            ## next item
        ## end def __iter__

        def __contains__(self, value):
            return (value in self.__dict__ and value == self.__dict__[value])
        ## end def __contains__
    ## end class __metaclass__

## end class OpanEnum


class EnumDispDirection(OpanEnum):
    """ Enumeration class for displacement directions.

    Contains enumeration parameters to indicate the displacement of the
    molecular geometry associated with a gradient, hessian or other object.

    **Enum Values**

    """

    # CHANGE OPANERROR METACLASS EXAMPLE IF ANY NEW TYPECODES ARE ADDED
    # HERE!!! (Should never happen, but....)

    #: Positive displacement along a particular normal mode
    POSITIVE = 'POSITIVE'

    #: Non-displaced geometry
    NO_DISP = 'NO_DISP'

    #: Negative displacement along a particular normal mode
    NEGATIVE = 'NEGATIVE'

## end class EnumDispDirection


class EnumMassPertType(OpanEnum):
    """ Enumeration class for atom mass perturbation types.

    Contains enumeration parameters to indicate the type of mass perturbation
    to be applied to the various atoms of a geometry, in order to: (a) break
    inertial degeneracy sufficiently to allow VPT2 computation using a lower-
    symmetry formalism; and (b), in the case of linear molecules, introduce
    an artificial 'directional preference' to the masses of the atoms
    (:attr:`BY_COORD` enum value) to break the intrinsic degeneracy
    of the bending modes.

    **Enum Values**

    """

    #: Atomic masses are used without modification
    NO_PERTURB = 'NO_PERTURB'

    #: Atomic masses are perturbed atom-by-atom in an isotropic fashion
    BY_ATOM = 'BY_ATOM'

    #: Atomic masses are perturbed anisotropically, where the
    #: perturbation factor for each atom's mass varies slightly in the
    #: x-, y-, and z-directions
    BY_COORD = 'BY_COORD'

## end class EnumMassPertType


class EnumTopType(OpanEnum):
    """ Enumeration class for classifying types of molecular tops.

    Contains enumeration parameters to indicate the type of molecular top
    associated with a particular geometry.

    Inertial moments with magnitudes less than
    :attr:`opan.const.PRM.ZERO_MOMENT_TOL` are taken as zero. Nonzero moments
    by this metric are considered to be equal if their ratio differs from
    unity by less than :attr:`opan.const.PRM.EQUAL_MOMENT_TOL`.  See
    :func:`opan.utils.inertia.principals` and the other functions defined
    in :mod:`opan.utils.inertia` for more details.

    **Enum Values**

    """

    #: Three zero principal inertial moments
    ATOM = 'ATOM'

    #: One zero and two equal non-zero moments
    LINEAR = 'LINEAR'

    #: Three equal, non-zero moments
    SPHERICAL = 'SPHERICAL'

    #: Three non-zero moments; largest two equal
    SYMM_PROL = 'SYMM_PROL'

    #: Three non-zero moments; smallest two equal
    SYMM_OBL = 'SYMM_OBL'

    #: Three unique non-zero moments
    ASYMM = 'ASYMM'

## end class EnumTopType


class EnumSoftware(OpanEnum):
    """ Enumeration class for identifying computational chemistry packages.

    This enum will be expanded if/when support for additional packages is
    implemented.

    **Enum Values**

    """

    #: The |orca| program package
    ORCA = 'ORCA'

## end class EnumSoftware


class EnumFileType(OpanEnum):
    """ Enumeration class for the file types generated by computational codes.

    **Enum Values**

    """

    #: XYZ atomic coordinates, assumed to follow the `Open Babel XYZ
    #: specification <http://openbabel.org/docs/2.3.0/FileFormats/
    #: XYZ_cartesian_coordinates_format.html>`_ |external link|
    XYZ = 'XYZ'

    #: Files containing nuclear gradient information
    GRAD = 'GRAD'

    #: Files containing nuclear Hessian information
    HESS = 'HESS'

    #: Files containing computational output
    OUTPUT = 'OUTPUT'

    #: Input files for defining computations
    INPUTFILE = 'INPUTFILE'

## end class EnumFileType


class EnumAnharmRepoData(OpanEnum):
    """ Enumeration class for datatypes in VPT2 HDF5 repository.

    Contains enumeration parameters to indicate the type of data to be
    retrieved from the on-disk HDF5 repository in the VPT2 anharmonic
    calculations of the :mod:`opan.anharm` submodule.

    **Enum Values**

    """

    #: Energy value at the given displacement
    ENERGY = 'ENERGY'

    #: Geometry at the given displacement (max precision available)
    GEOM = 'GEOM'

    #: Cartesian gradient vector
    GRAD = 'GRAD'

    #: Cartesian Hessian matrix
    HESS = 'HESS'

## end class EnumAnharmRepoData


class EnumAnharmRepoParam(OpanEnum):
    """ Enumeration class for parameters in VPT2 HDF5 repository.

    Contains enumeration parameters to indicate the parameter value to be
    retrieved from the on-disk HDF5 repository in the VPT2 anharmonic
    calculations of the :mod:`opan.anharm` submodule.

    **Enum Values**

    """

    #: Length-`N` vector of all-caps atomic symbols
    ATOMS = 'ATOMS'

    #: Displacement increment in :math:`\mathrm{B}\,\mathrm{u^{1/2}}`.
    #: Note that these are *not* atomic |units|, which would instead be
    #: :math:`\mathrm{B}\,\mathrm{m_e^{1/2}}`.
    INCREMENT = 'INCREMENT'

    #: Cartesian center of mass of system, *with perturbations applied*,
    #: if any
    CTR_MASS = 'CTR_MASS'

    #: Reference values of atomic masses (unperturbed)
    REF_MASSES = 'REF_MASSES'

    #: :class:`opan.const.EnumMassPertType` indicating perturbation type
    PERT_MODE = 'PERT_MODE'

    #: Length-`3*N` vector of perturbation factors (all should be ~1)
    PERT_VEC = 'PERT_VEC'

## end class EnumAnharmRepoParam


class EnumUnitsRotConst(OpanEnum):
    """ Units Enumeration class for rotational constants.

    Contains enumeration parameters to indicate the associated/desired units
    of interpretation/display of a rotational constant.

    String expressions of these units are provided in
    :attr:`UNITS.rotConst`.

    #DOC: Add link to exposition(?) of how RotConst expression is developed,
    once written.

    **Enum Values**

    """

    #: Inverse moment of inertia, :math:`\frac{1}{\mathrm{u\,B^2}}`. Note that
    #: the mass |units| here are *not* atomic units, which would require
    #: :math:`\frac{1}{\mathrm{m_e\,B^2}}`.
    INV_INERTIA = 'INV_INERTIA'

    #: Angular frequency in atomic |units|, :math:`\frac{1}{\mathrm{T_a}}`
    ANGFREQ_ATOMIC = 'ANGFREQ_ATOMIC'

    #: Angular frequency in SI units, :math:`\frac{1}{\mathrm s}`
    #: (**NOT** :math:`\mathrm{Hz}`!)
    ANGFREQ_SECS = 'ANGFREQ_SECS'

    #: Cyclic frequency in atomic |units|,
    #: :math:`\frac{\mathrm{cyc}}{\mathrm{T_a}}`
    CYCFREQ_ATOMIC = 'CYCFREQ_ATOMIC'

    #: Cyclic frequency in :math:`\mathrm{Hz}`,
    #: :math:`\frac{\mathrm{cyc}}{\mathrm s}`
    CYCFREQ_HZ = 'CYCFREQ_HZ'

    #: Cyclic frequency in :math:`\mathrm{MHz}`,
    #: millions of :math:`\frac{\mathrm{cyc}}{\mathrm s}`
    CYCFREQ_MHZ = 'CYCFREQ_MHZ'

    #: Wavenumbers in atomic |units|, :math:`\frac{\mathrm{cyc}}{\mathrm{B}}`
    WAVENUM_ATOMIC = 'WAVENUM_ATOMIC'

    #: Wavenumbers in conventional units,
    #: :math:`\frac{\mathrm{cyc}}{\mathrm{cm}}`
    WAVENUM_CM = 'WAVENUM_CM'

## end class EnumUnitsRotConst


# ======  Constants Classes  ====== #

class CIC(object):
    """Container for application-internal code information constants

    These may need expansion into dictionaries keyed by
    :class:`~opan.const.EnumSoftware` enum values, depending on the atoms
    supported by various software packages. They may also require adjustment
    to accommodate 'special' atom types such as ghost atoms and point charges.

    **Members**

    """

    #: |int| -- Maximum atomic number supported
    MAX_ATOMIC_NUM = 103

    #: |int| -- Minimum atomic number supported
    MIN_ATOMIC_NUM = 1

## end class CIC


class PHYS(object):

    """Container for physical constants

    **Members**

    """

    # Imports
    import numpy as _np

    #: |float| --
    #: Angstroms per Bohr radius (source: `NIST <http://physics.nist.gov/
    #: cgi-bin/cuu/Value?bohrrada0|search_for=bohr+radius>`__ |external link|)
    ANG_PER_BOHR = 0.52917721067

    #: |float| --
    #: Electron mass per unified atomic mass unit (source: `NIST
    #: <http://physics.nist.gov/cgi-bin/cuu/Value?meu|
    #: search_for=electron+mass>`__ |external link|)
    ME_PER_AMU = 1822.8885

    #: |float| --
    #: Seconds per atomic time unit (source: `NIST <http://physics.nist.gov/
    #: cgi-bin/cuu/Value?aut|search_for=atomic+time+unit>`__ |external link|)
    SEC_PER_TA = 2.4188843265e-17

    #: |float| --
    #: Speed of light in atomic |units|, :math:`\frac{B}{T_a}`. Calculated from
    #: the `NIST <http://physics.nist.gov/cgi-bin/cuu/Value
    #: ?c|search_for=speed+of+light>`__ |external link| value for the speed of
    #: light in vacuum, :math:`2.99792458e8\ \frac{m}{s}`, using
    #: :attr:`ANG_PER_BOHR` and :attr:`SEC_PER_TA` as conversion factors
    LIGHT_SPEED = 137.036

    #: |float| --
    #: Standard Planck constant, equal to :math:`2\pi` in atomic |units| of
    #: :math:`\frac{\mathrm{E_h\,T_a}}{\mathrm{cyc}}`
    PLANCK = 2 * _np.pi

    #: |float| --
    #: Reduced Planck constant, unity by definition in the atomic |units|
    #: of :math:`\mathrm{E_h\,T_a}`
    PLANCK_BAR = 1

## end class PHYS


class DEF(object):
    """Container for default parameter values (possibly user-adjustable)
    """

    from .const import EnumSoftware as _E_SW, EnumFileType as _E_FT

    #: |float| --
    #: Relative magnitude of atomic mass perturbations
    MASS_PERT_MAG = 1e-4

    #Moment_Tol = 1e-4

    #: |float| --
    #: Acceptable deviation from Kronecker delta for orthonormality testing
    ORTHONORM_TOL = 1e-8

    #: |float| --
    #: Max precision of HESS geometries (currently |orca|-specific)
    HESS_COORD_MATCH_TOL = 1e-6

    #: |float| --
    #: Max precision of freqs in IR spectrum block (currently |orca|-specific)
    HESS_IR_MATCH_TOL = 1e-2

    #: |float| --
    #: Max precision of GRAD geometries (currently |orca|-specific)
    GRAD_COORD_MATCH_TOL = 1e-7

    #: |float| --
    #: Max tolerable deviation between XYZ geoms (currently |orca|-specific)
    XYZ_COORD_MATCH_TOL = 1e-12

    # |float| --
    # Required quality of coordinate match for symmetry detection
    SYMM_MATCH_TOL = 1e-3

    # |float| --
    # Tolerance for deviation in searching for neighbor axes in cubic
    # symmetry groups -- value is in **radians**, and equals 1.5 degrees
    SYMM_AXIS_MATCH_TOL = 0.026179939

    # |int| --
    # Rounding atomic masses to avoid precision errors in atom matching
    SYMM_ATWT_ROUND_DIGITS = 4

    # |int| --
    # Initial order of rotational symmetry to test for (conservative)
    SYMM_MATCH_NMAX = 10

    # |int| --
    # Maximum order of atom averaging when looking for possible axes of
    # rotational symmetry
    SYMM_AVG_MAX = 2

    #: |dict| of |dict| --
    #: Dictionary of dictionaries of file extensions for geometry, gradient,
    #: and hessian files from the various softwares.
    #:
    #: Access as :samp:`FILE_EXTS[{EnumSoftware}][{EnumFileType}]`
    FILE_EXTS = {
            _E_SW.ORCA :
                { _E_FT.GRAD : 'engrad',
                    _E_FT.HESS : 'hess',
                    _E_FT.XYZ : 'XYZ',
                    _E_FT.OUTPUT : 'out',
                    _E_FT.INPUTFILE : 'txt'
                    }
            }

## end class DEF


class PRM(object):
    """Container for internal computation parameters (not user-adjustable)

    **Members**

    """

    #: |float| --
    #: Minimum angle deviation (degrees) required for two vectors to be
    #: considered non-parallel
    NON_PARALLEL_TOL = 1e-3

    #: |float| --
    #: Vector magnitude below which a vector is considered equal to the zero
    #: vector; dimensionless or |units| of :math:`\mathrm{B}`
    ZERO_VEC_TOL = 1e-6

    #: |float| --
    #: Trap value for aberrantly large dipole  derivative values in |orca|
    #: if dipoles are not calculated in a NUMFREQ run
    MAX_SANE_DIPDER = 100.0

    #: |float| --
    #: Minimum deviation-ratio from unity below which two principal inertial
    #: moments are considered equal
    EQUAL_MOMENT_TOL = 1e-3

    #: |float| --
    #: Threshold value below which moments are  considered equal to zero;
    #: |units| of :math:`\mathrm{u\,B^2}`
    ZERO_MOMENT_TOL = 1e-3

## end class PRM


class SYMM(object):
    1 # Dummy line as first line to detach the docstring
    # Be sure to also re-doc-comment the parameters in DEF and de-comment the
    # SYMM class reference in the module docstring, in addition to completing
    # this docstring, once the symmetry tools are in workable shape.
    """ Container for constants used in symmetry determination.

    Angles are in radians.

    *UNDER DEVELOPMENT*

    TD_C2_C2 = OH_C4_C4 = 1.570796327 (90 deg)
    TD_C2_C3 = OH_C4_C3 = 0.955316618 (54.736 deg)
    TD_C3_C3 = OH_C3_C3 = 1.230959417 (70.529 deg)

    Icosahedral nomenclature is IH_Cn_xCm
        n,m are orders of principal axes
        x indicates the x-th nearest Cm axis to any given Cn axis
    IH_C5_1C5 = 1.107147967 (63.435 deg)
    IH_C5_2C5 = 2.034442786 (116.565 deg)

    IH_C5_1C3 = 0.651981429 (37.356 deg)
    IH_C5_2C3 = 1.381327787 (79.144 deg)
    IH_C5_3C3 = 1.760264867 (100.856 deg)

    IH_C5_1C2 = 0.557442890 (31.939 deg)
    IH_C5_2C2 = 1.024395943 (58.694 deg)
    IH_C5_3C2 = 1.515997158 (86.860 deg)

    IH_C3_1C3 = 0.729350896 (41.789 deg)
    IH_C3_2C3 = 1.230262025 (70.489 deg)
    IH_C3_3C3 = 1.910129997 (109.442 deg)

    IH_C3_1C2 = 0.364698782 (20.896 deg)
    IH_C3_2C2 = 0.957758721 (54.876 deg)
    IH_C3_3C2 = 1.592672312 (91.253 deg)

    IH_C2_1C2 = 0.632381404 (36.233 deg)
    IH_C2_2C2 = 1.054424138 (60.414 deg)
    IH_C2_3C2 = 1.264759664 (72.465 deg)
    IH_C2_4C2 = 1.467100515 (84.059 deg)

    """

    TD_C2_C2 = OH_C4_C4 = 1.570796327
    TD_C2_C3 = OH_C4_C3 = 0.955316618
    TD_C3_C3 = OH_C3_C3 = 1.230959417

    IH_C5_1C5 = 1.107147967
    IH_C5_2C5 = 2.034442786

    IH_C5_1C3 = 0.651981429
    IH_C5_2C3 = 1.381327787
    IH_C5_3C3 = 1.760264867

    IH_C5_1C2 = 0.557442890
    IH_C5_2C2 = 1.024395943
    IH_C5_3C2 = 1.515997158

    IH_C3_1C3 = 0.729350896
    IH_C3_2C3 = 1.230262025
    IH_C3_3C3 = 1.910129997

    IH_C3_1C2 = 0.364698782
    IH_C3_2C2 = 0.957758721
    IH_C3_3C2 = 1.592672312

    IH_C2_1C2 = 0.632381404
    IH_C2_2C2 = 1.054424138
    IH_C2_3C2 = 1.264759664
    IH_C2_4C2 = 1.467100515

## end class SYMM


class UNINIT(object):
    """Container for numerical values indicating an un-initialized parent
        object.

    #TODO: (?)UNINIT: Consider deprecating in favor of a custom Exception.

    **Members**
    """

    # Empty doc comments trigger inclusion in the documentation
    #: |int| --
    UNSIGNED_INT = -1

    #: |float| --
    UNSIGNED_FLOAT = -1.0

    #: |int| --
    SIGNED_INT = 1234567890

    #: |float| --
    SIGNED_FLOAT = 1234567890.12345

## end class UNINIT


class UNITS(object):
    """Container for dicts providing strings describing the various display
    units available for physical quantities.

    Dictionary keys are the enum values provided in the corresponding
    ``EU_xxxx`` class in this module (:mod:`opan.const`).

    ================== =========================== =====================
      Dictionary              Enum                  Physical Quantity
    ================== =========================== =====================
     :attr:`rotConst`   :attr:`EnumUnitsRotConst`   Rotational constant
    ================== =========================== =====================

    """

    #TODO: (as occurs) const.UNITS: Add dicts for other units-based enums
    #TODO: (?) If good way of prettyprinting the dict becomes available,
    #   add doc comments.

    # Imports
    from .const import EnumUnitsRotConst as _EUrc

    #: |dict| --
    rotConst = {
            _EUrc.INV_INERTIA :        "1/(amu*B^2)",
            _EUrc.ANGFREQ_ATOMIC :     "1/Ta",
            _EUrc.ANGFREQ_SECS :    "1/s",
            _EUrc.CYCFREQ_ATOMIC :  "cyc/Ta",
            _EUrc.CYCFREQ_HZ :      "cyc/s",
            _EUrc.CYCFREQ_MHZ :     "MHz",
            _EUrc.WAVENUM_ATOMIC :     "cyc/B",
            _EUrc.WAVENUM_CM :         "cyc/cm"
                }

## end class UNITS


#TODO:(?) const: do atom_num/Sym need to handle ghost atoms / charges?
# Must use str.upper() when retrieving values from atom_num dictionary
atom_num = {"H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
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

atom_sym = {1: 'H', 2: 'HE', 3: 'LI', 4: 'BE', 5: 'B', 6: 'C', 7: 'N',
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
