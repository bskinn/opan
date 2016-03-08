#-------------------------------------------------------------------------------
# Name:        opan_utils_inertia
# Purpose:     Test objects for opan.utils.inertia
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     7 Mar 2016
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


import unittest


class SuperOpanUtilsInertia(object):
    # Superclass for all test cases for the various top types. Implemented
    #  as a top-level class rather than a subclass of unittest.TestCase to
    #  prevent the tests being run spuriously on the superclass. Also avoids
    #  run-time errors since self.fname does not exist in the superclass
    #  namespace.
    #
    # Subclasses must also subclass unittest.TestCase, and must declare
    #  SuperOpanUtilsInertia as a higher priority superclass. Otherwise, the
    #  default form of the `setUpClass` method from unittest.TestCase will be
    #  called.

    # Imports
    import os

    # Common constants
    filedir = os.path.join('resource','inertia')

    # Common setUpClass method
    @classmethod
    def setUpClass(cls):

        # Imports
        import os
        from opan.xyz import OpanXYZ as XYZ
        from opan.hess import OrcaHess as HESS

        # Attach class objects
        basename = os.path.join(cls.filedir, cls.fname)
        cls.xyz = XYZ(path=basename + '.xyz')
        cls.hess = HESS(path=basename + '.hess')

        # Always long messages
        cls.longMessage = True

    def test_ctr_mass(self):
        import opan.utils.inertia as oui
        ctr_mass = oui.ctr_mass(self.xyz.geoms[0], self.hess.atom_masses)
        for i in range(3):
            self.assertAlmostEqual(self.ctr_mass[i], ctr_mass[i],
                        delta=1e-7,
                        msg="Center-of-mass index '" + str(i) + "'")

    def test_ctr_geom(self):
        import opan.utils.inertia as oui
        ctr_geom = oui.ctr_geom(self.xyz.geoms[0], self.hess.atom_masses)
        for i in range(ctr_geom.shape[0]):
            self.assertAlmostEqual(self.ctr_geom[i], ctr_geom[i],
                        delta=1e-7,
                        msg="Centered geometry index '" + str(i) + "'")

    def test_i_tensor(self):
        import opan.utils.inertia as oui
        i_tensor = oui.inertia_tensor(self.xyz.geoms[0],
                                                    self.hess.atom_masses)
        for i in range(i_tensor.shape[0]):
            for j in range(i_tensor.shape[1]):
                if i_tensor[i,j] < 1e-10:
                    self.assertAlmostEqual(self.i_tensor[i,j],
                            i_tensor[i,j],
                            delta=1e-7,
                            msg="Inertia tensor element (" + str(i) + "," +
                                                             str(j) + ")")
                else:
                    self.assertAlmostEqual(self.i_tensor[i,j] / i_tensor[i,j],
                            1.0, delta=1e-7,
                            msg="Inertia tensor element (" + str(i) + "," +
                                                             str(j) + ")")

    def test_moments(self):
        import opan.utils.inertia as oui
        moments = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[0]
        for i in range(moments.shape[0]):
            if moments[i] < 1e-10:
                self.assertAlmostEqual(self.moments[i], moments[i],
                        delta=1e-7,
                        msg="Principal moment index '" + str(i) + "'")
            else:
                self.assertAlmostEqual(self.moments[i] / moments[i], 1.0,
                        delta=1e-7,
                        msg="Principal moment index '" + str(i) + "'")

    def test_axes(self):
        import opan.utils.inertia as oui
        axes = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[1]
        for i in range(axes.shape[0]):
            for j in range(axes.shape[1]):
                if abs(axes[i,j]) < 1e-10:
                    self.assertAlmostEqual(self.axes[i,j], axes[i,j],
                            delta=1e-7,
                            msg="Principal axis #" + str(j) + ", element " +
                                                                        str(i))
                else:
                    self.assertAlmostEqual(self.axes[i,j] / axes[i,j], 1.0,
                            delta=1e-6,
                            msg="Principal axis #" + str(j) + ", element " +
                                                                        str(i))

    def test_axes_orthonorm(self):
        import opan.utils.inertia as oui
        from opan.utils.vector import orthonorm_check as onchk
        axes = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[1]
        on, nfail, ofail = onchk(axes, report=True)
        self.assertTrue(on, msg="Norm failures: " + str(nfail) +
                                "; ortho failures: " + str(ofail))

    def test_rot_consts(self):
        import opan.utils.inertia as oui
        from opan.const import PRM, EnumUnitsRotConst as EURC
        rc = oui.rot_consts(self.xyz.geoms[0], self.hess.atom_masses,
                                                            EURC.INV_INERTIA)
        for i in range(rc.shape[0]):
            if rc[i] >= 1/(2.0*PRM.ZERO_MOMENT_TOL):
                self.assertAlmostEqual(self.rc[i], rc[i],
                        delta=1.0,
                        msg="Rotational constant index '" + str(i) + "'")
            else:
                self.assertAlmostEqual(self.rc[i] / rc[i], 1.0,
                        delta=PRM.EQUAL_MOMENT_TOL,
                        msg="Rotational constant index '" + str(i) + "'")

    def test_toptype(self):
        import opan.utils.inertia as oui
        top = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[2]
        self.assertEqual(top, self.top)

## end class SuperOpanUtilsInertia


class TestOpanUtilsInertiaAsymm(SuperOpanUtilsInertia, unittest.TestCase):
    # Asymmetric molecule test-case (H2O). Includes proofreading tests for,
    #  e.g., bad object shapes

    # Imports
    import numpy as np
    from opan.const import EnumTopType, EnumUnitsRotConst as EURC

    # Constants for superclass method use
    fname = "H2O_Asymm"
    ctr_mass = np.array([-11.42671795,   0.2152039 ,  -0.02145012])
    ctr_geom = np.array([-0.07503316, -0.09586698,  0.04549049,
                1.7685783 ,  0.01137086, -0.00539566,
                -0.57765012,  1.51023212, -0.71663048])
    i_tensor = np.array([[ 2.99702139,  0.74400967, -0.35304507],
                [ 0.74400967,  4.13012039,  1.16077065],
                [-0.35304507,  1.16077065,  6.02553148]])
    moments = np.array([ 2.41469379,  4.16164284,  6.57633663])
    axes = np.array([[  8.16492638e-01,  -5.77355846e-01,   6.07765301e-15],
                [ -5.21610191e-01,  -7.37657520e-01,  -4.28700585e-01],
                [  2.47512789e-01,   3.50030872e-01,  -9.03446627e-01]])
    top = EnumTopType.ASYMM
    rc = np.array([ 0.207065592447,  0.120144860861,  0.076030171231])

    rc_units = {EURC.ANGFREQ_ATOMIC: np.array([  1.135920230159e-04,   6.590905634729e-05,   4.170862410476e-05]),
         EURC.ANGFREQ_SECS: np.array([  4.696050231566e+12,   2.724770904719e+12,   1.724291800473e+12]),
         EURC.CYCFREQ_ATOMIC: np.array([  1.807873195879e-05,   1.048975211219e-05,   6.638133695834e-06]),
         EURC.CYCFREQ_HZ: np.array([  7.473996073616e+11,   4.336607582790e+11,   2.744295633781e+11]),
         EURC.CYCFREQ_MHZ: np.array([ 747399.60736156872 ,  433660.758279046626,  274429.563378054823]),
         EURC.INV_INERTIA: np.array([ 0.207065592447,  0.120144860861,  0.076030171231]),
         EURC.WAVENUM_ATOMIC: np.array([  1.318826101078e-07,   7.652173233680e-08,   4.842454659134e-08]),
         EURC.WAVENUM_CM: np.array([ 24.922201369646,  14.460511669382,   9.150913076387])}


    def test_UtilsInertiaBadGeomShape(self):
        from opan.utils.inertia import ctr_mass as cm, _fadn_orth as fO, \
                                                        _fadn_par as fP
        import numpy as np
        self.assertRaises(ValueError, cm, np.array([[1,1],[1,1]]),
                                                    self.hess.atom_masses)
        self.assertRaises(ValueError, fO, [1,2,3], self.xyz.geoms[0][:-2])
        self.assertRaises(ValueError, fP, [1,2,3], self.xyz.geoms[0][:-2])

    def test_UtilsInertiaBadMassesShape(self):
        from opan.utils.inertia import ctr_mass as cm
        import numpy as np
        self.assertRaises(ValueError, cm, self.xyz.geoms[0],
                                                     np.array([[1,1],[1,1]]))

    def test_UtilsInertiaBadGeomLength(self):
        from opan.utils.inertia import ctr_mass as cm
        import numpy as np
        self.assertRaises(ValueError, cm, self.xyz.geoms[0][:-2],
                                                    self.hess.atom_masses)

    def test_UtilsInertiaBadGeomMassesLengths(self):
        from opan.utils.inertia import ctr_mass as cm
        import numpy as np
        self.assertRaises(ValueError, cm, self.xyz.geoms[0],
                                                self.hess.atom_masses[:-2])

    def test_UtilsInertiaBadRefVecShape(self):
        from opan.utils.inertia import _fadn_orth as fO, _fadn_par as fP
        import numpy as np
        self.assertRaises(ValueError, fO, [1,2,3,4], self.xyz.geoms[0])
        self.assertRaises(ValueError, fP, [1,2,3,4], self.xyz.geoms[0])

    def test_UtilsInertiaBadRefVecNorm(self):
        from opan.utils.inertia import _fadn_orth as fO, _fadn_par as fP
        import numpy as np
        self.assertRaises(ValueError, fO, [1e-7, 0, 1e-9], self.xyz.geoms[0])
        self.assertRaises(ValueError, fP, [1e-7, 0, 1e-9], self.xyz.geoms[0])

    def test_UtilsInertiaCheckAllRotConstUnits(self):
        from opan.utils.inertia import rot_consts
        from opan.const import EnumUnitsRotConst as EURC
        for (u, a) in [(u, rot_consts(self.xyz.geoms[0],
                    self.hess.atom_masses, units=u)) for u in EURC]:
            for i in range(a.shape[0]):
                self.assertAlmostEqual(self.rc_units[u][i] / a[i], 1.0,
                        delta=1e-3,
                        msg="Rotational constant units '" + str(u) +
                                                ",' index '" + str(i) + "'")

    def test_UtilsInertiaBadRotConstUnits(self):
        from opan.utils.inertia import rot_consts
        self.assertRaises(ValueError, rot_consts, self.xyz.geoms[0],
                                self.hess.atom_masses, units="ThisIsInvalid")

## end class TestOpanUtilsInertiaAsymm


class TestOpanUtilsInertiaAtom(SuperOpanUtilsInertia, unittest.TestCase):
    # Lone atom test-case (Cu).

    # Imports
    import numpy as np
    from opan.const import EnumTopType, PRM

    # Constants for superclass method use
    fname = "Cu_Atom"
    ctr_mass = np.array([-1.88972612,  3.77945225, -7.5589045 ])
    ctr_geom = np.array([0.0, 0.0, 0.0])
    i_tensor = np.zeros((3,3))
    moments = np.zeros((3,))
    axes = np.eye(3)
    top = EnumTopType.ATOM
    rc = np.repeat(1.0/(2.0*PRM.ZERO_MOMENT_TOL), 3)

## end class TestOpanUtilsInertiaAtom


class TestOpanUtilsInertiaLinear(SuperOpanUtilsInertia, unittest.TestCase):
    # Linear molecule test-case (chloroethyne).

    # Imports
    import numpy as np
    from opan.const import EnumTopType, PRM
    from opan.tests import assertErrorAndTypecode

    # Constants for superclass method use
    fname = "HC2Cl_Linear"
    ctr_mass = np.array([ 3.78157026, -0.13914054,  0.11918503])
    ctr_geom = np.array([ -1.81226617e+00,   4.52442528e-10,   7.10376657e-09,
                         1.29696963e+00,   4.73253989e-09,  -1.02091483e-08,
                         3.58244725e+00,  -6.42275907e-09,  -1.13833863e-08,
                         5.59880599e+00,   4.22696256e-09,   7.43858825e-09])
    i_tensor = np.array([[  5.44288079e-15,   2.07854690e-07,   1.06328882e-06],
                       [  2.07854690e-07,   3.22388376e+02,  -4.43485080e-16],
                       [  1.06328882e-06,  -4.43485080e-16,   3.22388376e+02]])
    moments = np.array([  5.68434189e-14,   3.22388376e+02,   3.22388376e+02])
    axes = np.array([[ -1.00000000e+00,   2.49655671e-10,  -3.91982518e-09],
                   [  2.49655671e-10,   1.00000000e+00,   0.00000000e+00],
                   [  3.91982518e-09,  -9.78606588e-19,  -1.00000000e+00]])
    top = EnumTopType.LINEAR
    rc = np.array([1.0/(2.0*PRM.EQUAL_MOMENT_TOL),
                                    1.550924404326e-03, 1.550924404326e-03])

    def test_UtilsInertiaLinearNoNonParallelVec(self):
        from opan.utils.inertia import _fadn_par as fP, ctr_geom as cg
        import numpy as np
        from opan.error import InertiaError as INErr
        self.assertErrorAndTypecode(INErr, fP , INErr.BAD_GEOM,
                            self.xyz.displ_single(0,0,1),
                            cg(self.xyz.geoms[0], self.hess.atom_masses))

## end class TestOpanUtilsInertiaLinear


class TestOpanUtilsInertiaSymmProl(SuperOpanUtilsInertia, unittest.TestCase):
    # Prolate symmetric test case (chloromethane)

    # Imports
    import numpy as np
    from opan.const import EnumTopType

    # Constants for superclass method use
    fname = "CH3Cl_SymmProl"
    ctr_mass = np.array([ 1.19526288,  1.19526288,  1.19526288])
    ctr_geom = np.array([-1.36880665, -1.36880665, -1.36880665,
                            0.61213832,  0.61213832, 0.61213832,
                            -2.54039904, -2.54039904, -0.13884799,
                            -0.13884799, -2.54039904, -2.54039904,
                            -2.54039904, -0.13884799, -2.54039905])
    i_tensor = np.array([[ 97.63769388, -43.0052599 , -43.00525991],
                       [-43.0052599 ,  97.63769391, -43.00525989],
                       [-43.00525991, -43.00525989,  97.63769389]])
    moments = np.array([  11.6271741 ,  140.64295379,  140.6429538 ])
    axes = np.array([[ -5.77350269e-01,  -4.08248290e-01,  -7.07106781e-01],
                       [ -5.77350269e-01,  -4.08248291e-01,   7.07106781e-01],
                       [ -5.77350269e-01,   8.16496581e-01,   2.78286699e-10]])
    top = EnumTopType.SYMM_PROL
    rc = np.array([ 0.043002710352,  0.003555101671,  0.003555101671])

## end class TestOpanUtilsInertiaSymmProl


class TestOpanUtilsInertiaSymmObl(SuperOpanUtilsInertia, unittest.TestCase):
    # Oblate symmetric test case (ammonia)

    # Imports
    import numpy as np
    from opan.const import EnumTopType

    # Constants for superclass method use
    fname = "NH3_SymmObl"
    ctr_mass = np.array([ -1.36098776e+01,   1.70206822e+00,  -7.62922335e-04])
    ctr_geom = np.array([ -4.79117194e-02,  -1.35356483e-01,   5.88277405e-03,
                         1.89623378e+00,   3.55505121e-02,  -1.54206207e-03,
                        -6.15221552e-01,   8.55891028e-01,  -1.57668545e+00,
                        -6.15238957e-01,   9.89449583e-01,   1.49648146e+00])
    i_tensor = np.array([[ 6.74683794,  0.98560398, -0.04281628],
                           [ 0.98560398,  9.1833784 , -0.12106317],
                           [-0.04281628, -0.12106317,  6.40284824]])
    moments = np.array([ 6.39758672,  6.39807449,  9.53740337])
    axes = np.array([[-0.941357532609, -0.052449650272, -0.333309210773],
                   [ 0.334894464181, -0.024822665809, -0.94192862422 ],
                   [ 0.041130203771, -0.998315015137,  0.040932100963]])
    top = EnumTopType.SYMM_OBL
    rc = np.array([ 0.078154470132,  0.078148511783,  0.05242517071 ])

## end class TestOpanUtilsInertiaSymmObl


class TestOpanUtilsInertiaPlanar(SuperOpanUtilsInertia, unittest.TestCase):
    # Planar (oblate symmetric special case) test case (benzene)

    # Imports
    import numpy as np
    from opan.const import EnumTopType

    # Constants for superclass method use
    fname = "C6H6_Planar"
    ctr_mass = np.array(
            [0.0, 0.0, 0.0])
    ctr_geom = np.array([ 1.889726124565,  0.            ,  0.,
                0.944863062283, 1.636550830098,  0.,
                -0.944863062283,  1.636550830098, 0.,
                -1.889726124565,  0.            ,  0.,
                -0.944863062283, -1.636550830098,  0.,
                0.944863062283, -1.636550830098,  0.,
                3.77945224913 ,  0., 0.,
                1.889726124565,  3.273101660137,  0.,
                -1.889726124565,  3.273101660137,  0.,
                -3.77945224913 , 0.,  0.,
                -1.889726124565, -3.273101660137, 0.,
                 1.889726124565, -3.273101660137,  0.])
    i_tensor = np.array([[ 171.871779008624,   -0.,   -0.],
                            [  -0.,  171.871779003993,   -0.],
                            [  -0., -0.,  343.743558012617]])
    moments = np.array([171.871779003993, 171.871779008624, 343.743558012617])
    axes = np.array([[ 1.,  0.,  0.],
                       [ 0.,  1.,  0.],
                       [ 0.,  0.,  1.]])
    top = EnumTopType.SYMM_OBL
    rc = np.array([ 0.002909145427,  0.002909145427,  0.001454572714])

## end class TestOpanUtilsInertiaPlanar


class TestOpanUtilsInertiaSpher(SuperOpanUtilsInertia, unittest.TestCase):
    # Spherical test case (methane)

    # Imports
    import numpy as np
    from opan.const import EnumTopType

    # Constants for superclass method use
    fname = "CH4_Spher"
    ctr_mass = np.array(
            [  8.734135537996e-11,   4.387156003093e-11,  -1.188484647673e-11])
    ctr_geom = np.array(
            [  4.000728815448e-11,   2.009566928560e-11,  -5.443942085537e-12,
             1.180503717314e+00,   1.180503717408e+00,   1.180503717381e+00,
            -1.180503717657e+00,  -1.180503717520e+00,   1.180503717479e+00,
             1.180503717435e+00,  -1.180503717544e+00,  -1.180503717454e+00,
            -1.180503717568e+00,   1.180503717417e+00,  -1.180503717341e+00])
    i_tensor = np.array(
            [[  1.123790191290e+01,   7.775846633251e-11,   5.015445836420e-10],
           [  7.775846633251e-11,   1.123790191310e+01,  -3.777778090353e-11],
           [  5.015445836420e-10,  -3.777778090353e-11,   1.123790191366e+01]])
    moments = np.array([ 11.237901912633,  11.237901913117,  11.237901913909])
    axes = np.array(
            [[  5.773502691632e-01,  -4.082482905123e-01,   7.071067811802e-01],
           [  5.773502692095e-01,  -4.082482904247e-01,  -7.071067811929e-01],
           [  5.773502691962e-01,   8.164965809231e-01,   6.946399011554e-11]])
    top = EnumTopType.SPHERICAL
    rc = np.array([ 0.044492290811,  0.044492290809,  0.044492290806])

## end class TestOpanUtilsInertiaSpher

def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanUtilsInertiaAsymm),
                tl.loadTestsFromTestCase(TestOpanUtilsInertiaAtom),
                tl.loadTestsFromTestCase(TestOpanUtilsInertiaLinear),
                tl.loadTestsFromTestCase(TestOpanUtilsInertiaPlanar),
                tl.loadTestsFromTestCase(TestOpanUtilsInertiaSpher),
                tl.loadTestsFromTestCase(TestOpanUtilsInertiaSymmObl),
                tl.loadTestsFromTestCase(TestOpanUtilsInertiaSymmProl)
                ])
    return s

## end def suite




if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

