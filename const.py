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


"""Module with container classes and dictionaries bearing assorted constants
    for OpenAnharmonic.


Constants Classes
-----------------
CIC     -- Application-internal code information constants
PHYS    -- Physical constants
DEF     -- Default values for (theoretically) user-adjustable parameters
PRM     -- Internal computation parameters, intended to be non-user-
             adjustable
UNINIT  -- Constants representing un-initialized values
UNITS   -- Functions returning text strings of units descriptions

Enumeration Classes
-------------------
E_DispDirection     -- Displacement direction along a particular mode
E_MassPerturbation  -- Type of atomic mass perturbation being applied
E_TopType           -- Molecular top classification

Units Enumeration Classes
-------------------------
EU_RotConst         -- Rotational constants

Methods
-------
keydef  -- Helper function for enumeration class variable definition

Dictionaries
------------
atomNum  -- Atomic number lookup from element symbol
                Note: keys are *all uppercase* (e.g., AR, NE, and NA)
atomSym  -- Element symbol lookup from atomic number
                Note: symbols are returned as *all uppercase*

Variables
---------
infty    -- Infinity symbol as Unicode string

"""


# Infinity symbol
infty = u"\u221E"


class E_DispDirection(object):
    """ Enumeration class for displacement directions.

    Contains enumeration parameters to indicate the displacement of the
    molecular geometry associated with a gradient, hess or other object

    Values are collected into the tuple 'E' for convenient Enum membership
    testing.

    Enum Values
    -----------
    Uninitialized : Object created without declaration of a displacement
            direction. [removed]
    Positive : Positive displacement along a particular normal mode.
    NoDisp : Non-displaced geometry.
    Negative : Negative displacement along a particular normal mode.
    """

    Positive = 'Positive'
    NoDisp = 'NoDisp'
    Negative = 'Negative'

    E = frozenset([
        Positive,
        NoDisp,
        Negative
        ])

## end class E_DispDirection


class E_MassPertType(object):
    """ Enumeration class for atom mass perturbation types.

    Contains enumeration parameters to indicate the type of mass perturbation
    to be applied to the various atoms of a geometry, in order to: (a) break
    inertial degeneracy sufficiently to allow VPT2 computation using a lower-
    symmetry formalism; and (b), in the case of linear molecules, introduce
    an artificial 'directional preference' to the masses of the atoms (ByCoord
    enum value) to break the intrinsic degeneracy of the bending modes.

    Values are collected into the tuple 'E' for convenient Enum membership
    testing.

    Enum Values
    -----------
    Uninitialized : Object created without declaration of a mass perturbation
            type. [removed]
    NoPerturb : Atomic masses are used without modification
    ByAtom : Atomic masses are perturbed atom-by-atom in an isotropic fashion.
    ByCoord : Atomic masses are perturbed anisotropically, where the
        perturbation factor for each atom's mass varies slightly in the
        x-, y-, and z-directions.
    """

    NoPerturb = 'NoPerturb'
    ByAtom = 'ByAtom'
    ByCoord = 'ByCoord'

    E = frozenset([
        NoPerturb,
        ByAtom,
        ByCoord
        ])

## end class E_MassPertType


class E_TopType(object):
    """ Enumeration class for classifying types of molecular tops.

    Contains enumeration parameters to indicate the type of molecular top
    associated with a particular geometry.

    Values are collected into the tuple 'E' for convenient Enum membership
    testing.

    Enum Values
    -----------
    Uninitialized : Object created without declaration of a top type.
                [removed]
    Invalid : Inertia matrix malformed during computation.
                (realistically, should never be needed) [removed]
    Atom : Three zero principal inertial moments.
    Linear : One zero and two equal non-zero moments.
    Spherical : Three equal, non-zero moments.
    SymmetricalProlate : Three non-zero moments; largest two equal.
    SymmetricalOblate : Three non-zero moments; smallest two equal.
    Asymmetrical : Three unique non-zero moments.
    """

    Atom = 'Atom'
    Linear = 'Linear'
    Spherical = 'Spherical'
    SymmProlate = 'SymmProlate'
    SymmOblate = 'SymmOblate'
    Asymmetrical = 'Asymmetrical'

    E = frozenset([
        Atom,
        Linear,
        Spherical,
        SymmProlate,
        SymmOblate,
        Asymmetrical
        ])

## end class E_TopType


class E_Software(object):
    """ #DOC: Docstring for E_Software

    """

    ORCA = 'ORCA'

    E = frozenset([
        ORCA
        ])

## end class E_Software


class E_FileType(object):
    """ #DOC: Docstring for E_FileType
    """

    xyz = 'xyz'
    grad = 'grad'
    hess = 'hess'
    output = 'output'
    inputfile = 'inputfile'

    E = frozenset([
        xyz,
        grad,
        hess,
        output,
        inputfile
        ])

## end class E_FileType


class ER_Data(object):
    """ Enumeration class for datatypes in HDF5 repository.

    Contains enumeration parameters to indicate the type of data to be
    retrieved from the on-disk HDF5 repository.

    Values are collected into the tuple 'E' for convenient Enum membership
    testing.

    Enum Values
    -----------
    en      : Energy value at the given displacement
    geom    : Geometry at the given displacement from the XYZ file
    grad    : Cartesian gradient vector as retrieved from the ENGRAD file
    hess    : Cartesian Hessian matrix as retrieved from the HESS file
    """

    en = 'en'
    geom = 'geom'
    grad = 'grad'
    hess = 'hess'

    E = frozenset([
        en,
        geom,
        grad,
        hess
        ])

## end class ER_Data


class ER_Param(object):
    """ Enumeration class for parameters stored in the HDF5 repository.

    Contains enumeration parameters to indicate the parameter value to be
    retrieved from the on-disk HDF5 repository.

    Values are collected into the tuple 'E' for convenient Enum membership
    testing.

    Enum Values
    -----------
    atoms       : Vector of all-caps atomic symbols
    increment   : Displacement increment in Bohr*amu^(1/2)
    ctr_mass    : Cartesian center of mass of *MASS-UNPERTURBED* system
    ref_masses  : Reference values of atomic masses (unperturbed)
    pert_mode   : E_MassPertType enum indicating perturbation type
    pert_vec    : 3*N vector of perturbation factors (all ~1)
    """

    atoms = 'atoms'
    increment = 'increment'
    ctr_mass = 'ctr_mass'
    ref_masses = 'ref_masses'
    pert_mode = 'pert_mode'
    pert_vec = 'pert_vec'
    E = frozenset([
        atoms,
        increment,
        ctr_mass,
        ref_masses,
        pert_mode,
        pert_vec
        ])

## end class ER_Param


class EU_RotConst(object):
    """ Units Enumeration class for rotational constants.

    Contains enumeration parameters to indicate the associated/desired units
    of interpretation/display of a rotational constant.

    Values are collected into the tuple 'E' for convenient Enum membership
    testing.

    More detailed exposition of the form/nature of the various units is
    provided in const.UNITS.rotConst

    Enum Values
    -----------
    InvInertia : Inverse moment of inertia
    AngFreqAtomic : Angular frequency, a.u., 1/Ta
    AngFreqSeconds : Angular frequency, 1/s
    CyclicFreqAtomic : Cyclic frequency, a.u., cyc/Ta
    CyclicFreqHz : Cyclic Frequency, Hz (cyc/s)
    CyclicFreqMHz : Cyclic Frequency, MHz (1e6 cyc/s)
    WaveNumAtomic : Wavenumber, a.u., cyc/Bohr
    WaveNumCM : Wavenumber, cyc/cm
    """

    InvInertia = 'InvInertia'
    AngFreqAtomic = 'AngFreqAtomic'
    AngFreqSeconds = 'AngFreqSeconds'
    CyclicFreqAtomic = 'CyclicFreqAtomic'
    CyclicFreqHz = 'CyclicFreqHz'
    CyclicFreqMHz = 'CyclicFreqMHz'
    WaveNumAtomic = 'WaveNumAtomic'
    WaveNumCM = 'WaveNumCM'

    E = frozenset([
        InvInertia,
        AngFreqAtomic,
        AngFreqSeconds,
        CyclicFreqAtomic,
        CyclicFreqHz,
        CyclicFreqMHz,
        WaveNumAtomic,
        WaveNumCM
        ])

## end class EU_RotConst


class CIC(object):
    """Container for application-internal code information constants
    Invalid_Atom_Symbol = -1
    Unsupported_Atomic_Number = "INVALID"
    Max_Atomic_Num = 103
    Min_Atomic_Num = 1
    """

    Invalid_Atom_Symbol = -1
    Unsupported_Atomic_Number = "INVALID"
    Max_Atomic_Num = 103
    Min_Atomic_Num = 1

## end class CIC


class PHYS(object):

    """Container for physical constants

    Ang_per_Bohr = 0.52917721092
    amu_per_me = 1822.8885
    sec_per_Ta = 2.4188843265e-17
    light_speed = 137.082  (B/Ta units)
    Planck = 6.28318531  (value of 2*pi; Eh*Ta/cyc units)
    Planck_bar = 1  (Eh*Ta units)
    """

    Ang_per_Bohr = 0.52917721092
    amu_per_me = 1822.8885
    sec_per_Ta = 2.4188843265e-17
    light_speed = 137.082
    Planck = 6.28318531
    Planck_bar = 1


class DEF(object):
    """Container for default parameter values (possibly user-adjustable)

    Mass_Perturbation_Magnitude = 1e-4 (relative to each atom's default mass)
    Moment_Tol = 1e-4  (Equality test tolerance for inertial moments)
    Orthonorm_Tol = 1e-8  (Acceptable deviation from Kronecker delta for
                            orthonormality testing)
    HESS_Coord_Match_Tol = 1e-6 (Max precision of HESS geometries)
    HESS_IR_Match_Tol = 1e-2 (Max precision of freqs in IR spectrum block)
    GRAD_Coord_Match_Tol = 1e-7 (Max precision of GRAD geometries)
    XYZ_Coord_Match_Tol = 1e-12 (Max tolerable deviation between XYZ geoms)
    Symm_Match_Tol = 1e-3 (for quality of fit in symmetry matching)
    Symm_Axis_Match_Tol = 0.026179939 (radians, equals 1.5 deg; tolerance for
                    deviation in searching for neighbor axes in cubic
                    symmetry groups)
    Symm_AtWt_Round_Digits = 4 (to avoid precision errors in atom matching)
    Symm_Match_nMax = 10 (initial order of symmetry rotation to test for)
    Symm_Avg_Max = 2 (maximum order of atom averaging when looking for
                        possible axes of rotational symmetry)
    File_Extensions (dictionary of dictionaries of file extensions for
                geom, gradient, and hessian files from the various softwares;
                access as F_E[E_Software][E_FileType])
    """

    from .const import E_Software as E_SW, E_FileType as E_FT

    Mass_Perturbation_Magnitude = 1e-4
    Moment_Tol = 1e-4
    Orthonorm_Tol = 1e-8
    HESS_Coord_Match_Tol = 1e-6
    HESS_IR_Match_Tol = 1e-2
    GRAD_Coord_Match_Tol = 1e-7
    XYZ_Coord_Match_Tol = 1e-12
    Symm_Match_Tol = 1e-3
    Symm_Axis_Match_Tol = 0.026179939
    Symm_AtWt_Round_Digits = 4
    Symm_Match_nMax = 10
    Symm_Avg_Max = 2

    File_Extensions = {
            E_SW.ORCA :
                { E_FT.grad : 'engrad',
                    E_FT.hess : 'hess',
                    E_FT.xyz : 'xyz',
                    E_FT.output : 'out',
                    E_FT.inputfile : 'txt'
                    }
            }

## end class DEF


class PRM(object):
    """Container for internal computation parameters (not user-adjustable)

    Non_Parallel_Tol = 1e-3
    Zero_Vec_Tol = 1e-6
    Symm_Min_Origin_Dist = 0.5
    Max_Sane_DipDer = 100
    """

    Non_Parallel_Tol = 1e-3
    Zero_Vec_Tol = 1e-6
    Symm_Min_Origin_Dist = 0.5
    Max_Sane_DipDer = 100

## end class PRM


class SYMM(object):
    """ Container for constants used in symmetry determination.

    Angles are in radians.

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


class UNINIT(object):
    """Container for numerical values indicating an un-initialized parent
        object.

    #TODO: (?)UNINIT: Consider deprecating in favor of a custom Exception.

    Unsigned_Long = -1
    Unsigned_Double = -1
    Signed_Long = 1234567890
    Signed_Double = 1234567890.12345
    """

    Unsigned_Long = -1
    Unsigned_Double = -1
    Signed_Long = 1234567890
    Signed_Double = 1234567890.12345


class UNITS(object):
    """Container for dicts providing strings describing the various display
        units available for physical quantities.

    Dictionary keys are the enum values provided in the corresponding
        EU_xxxx class in this module (const).


      Dictionary  | ENUM                | Physical Quantity
    -------------------------------------------------------------------------
       rotConst   |  .RotConstUnits     |  Rotational constant
    -------------------------------------------------------------------------
    """

    #TODO: (as occurs) const.UNITS: Add dicts for other units-based enums

    # Imports
    from .const import EU_RotConst as EUrc

    rotConst = {
            EUrc.InvInertia :        "1/(amu*B^2)",
            EUrc.AngFreqAtomic :     "1/Ta",
            EUrc.AngFreqSeconds :    "1/s",
            EUrc.CyclicFreqAtomic :  "cyc/Ta",
            EUrc.CyclicFreqHz :      "cyc/s",
            EUrc.CyclicFreqMHz :     "MHz",
            EUrc.WaveNumAtomic :     "cyc/B",
            EUrc.WaveNumCM :         "cyc/cm"
                }


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


if __name__ == '__main__':
    print("Module not executable.")
