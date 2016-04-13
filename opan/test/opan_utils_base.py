#-------------------------------------------------------------------------------
# Name:        opan_utils_base
# Purpose:     Test objects for opan.utils.base
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     28 Feb 2016
# Copyright:   (c) Brian Skinn 2016
# License:     The MIT License; see "license.txt" for full license terms
#                   and contributor agreement.
#
#       This file is part of opan (Open Anharmonic), a system for automated
#       computation of anharmonic properties of molecular systems via wrapper
#       calls to computational/quantum chemical software packages.
#
#       http://www.github.com/bskinn/opan
#
#-------------------------------------------------------------------------------


import unittest


class TestOpanUtilsBaseMisc(unittest.TestCase):

    def test_Utils_PackTupsGoodPacking(self):
        from opan.utils import pack_tups
        tups = pack_tups(range(3), range(3,6), range(6,9))
        [[self.assertEqual(tups[i][j], 3*j + i) for i in range(3)]
                                                for j in range(3)]

    def test_Utils_PackTupsStrNoIter(self):
        from opan.utils import pack_tups
        tups = pack_tups("ab", range(2))
        self.assertEqual(tups[0][0], "ab")
        self.assertEqual(tups[1][0], "ab")

    def test_Utils_PackTupsErrIfDiffLens(self):
        from opan.utils import pack_tups
        self.assertRaises(ValueError, pack_tups, range(2), range(3))

    def test_Utils_PackTupsTestAllSingletons(self):
        from opan.utils import pack_tups
        tups = pack_tups(0,1,2,3,4)
        self.assertEqual(len(tups), 1)
        self.assertTupleEqual(tups[0], tuple(range(5)))

    def test_Utils_SafeCastNumpyArray(self):
        import numpy as np
        from opan.utils import safe_cast as scast
        a = np.array(range(5))
        self.assertRaises(TypeError, scast, a, np.float_)

    def test_Utils_SafeCastIntToFloat(self):
        from opan.utils import safe_cast as scast
        v = scast(1, float)
        self.assertIsInstance(v, float)

    def test_Utils_MakeTimeStampSecs(self):
        from opan.utils import make_timestamp as mt
        self.assertEqual(mt(5), "0h 0m 5s")

    def test_Utils_MakeTimeStampMins(self):
        from opan.utils import make_timestamp as mt
        self.assertEqual(mt(500), "0h 8m 20s")

    def test_Utils_MakeTimeStampHours(self):
        from opan.utils import make_timestamp as mt
        self.assertEqual(mt(20000), "5h 33m 20s")

    def test_Utils_MakeTimeStampLongHours(self):
        from opan.utils import make_timestamp as mt
        self.assertEqual(mt(200000), "55h 33m 20s")

    def test_Utils_DeltaFxn(self):
        from opan.utils import delta_fxn as dfx

        self.assertEqual(dfx(1,1), 1)
        self.assertEqual(dfx(1,2), 0)


class TestOpanUtilsBaseCheckGeom(unittest.TestCase):
    import numpy as np

    coords = np.array(range(12))
    atoms = ['H', 'O', 'O', 'H']

    def test_Utils_CheckGeom_GoodCheck(self):
        from opan.utils import check_geom as cg

        # check_geom returns a tuple; success or not is the first element
        self.assertTrue(cg(self.coords, self.atoms,
                           self.coords, self.atoms)[0])

    def test_Utils_CheckGeom_ArgsNotVectors(self):
        from opan.utils import check_geom as cg
        import numpy as np

        self.assertRaises(ValueError, cg,
                          self.coords.reshape((3,4)), self.atoms,
                          self.coords, self.atoms)
        self.assertRaises(ValueError, cg,
                          self.coords, np.array(self.atoms).reshape((2,2)),
                          self.coords, self.atoms)
        self.assertRaises(ValueError, cg,
                          self.coords, self.atoms,
                          self.coords.reshape((3, 4)), self.atoms)
        self.assertRaises(ValueError, cg,
                          self.coords, self.atoms,
                          self.coords, np.array(self.atoms).reshape((2, 2)))

    def test_Utils_CheckGeom_CoordAtomSizeMismatch(self):
        from opan.utils import check_geom as cg

        self.assertRaises(ValueError, cg,
                          self.coords, self.atoms,
                          self.coords[:-1], self.atoms)
        self.assertRaises(ValueError, cg,
                          self.coords[:-1], self.atoms,
                          self.coords, self.atoms)

    def test_Utils_CheckGeom_GeomSizeMismatch(self):
        from opan.utils import check_geom as cg

        tup = cg(self.coords, self.atoms, self.coords[:9], self.atoms[:3])

        self.assertFalse(tup[0])
        self.assertEqual(tup[1], 'coord_dim_mismatch')
        self.assertIsNone(tup[2])

    def test_Utils_CheckGeom_CoordMismatch(self):
        from opan.utils import check_geom as cg

        # Change one of the coordinates to achieve mismatch
        coord_mod = self.coords.copy()
        coord_mod[0] = -12

        # Store the call result and check its contents
        tup = cg(coord_mod, self.atoms, self.coords, self.atoms)
        self.assertFalse(tup[0])
        self.assertEqual(tup[1], 'coord_mismatch')
        self.assertFalse(tup[2][0])
        self.assertTrue(all(tup[2][1:]))

    def test_Utils_CheckGeom_AtomMismatch(self):
        from opan.utils import check_geom as cg

        # Change one of the atoms to achieve mismatch
        atoms_mod = self.atoms.copy()
        atoms_mod[2] = 'C'

        # Store the call result and check its contents
        tup = cg(self.coords, atoms_mod, self.coords, self.atoms)
        self.assertFalse(tup[0])
        self.assertEqual(tup[1], 'atom_mismatch')
        self.assertFalse(tup[2][2])
        self.assertTrue(all(tup[2][:2]))
        self.assertTrue(tup[2][3])


class TestOpanUtilsBaseTemplateSubst(unittest.TestCase):
    from textwrap import dedent

    template = dedent("""\
        This <NOUN> is a two-line string.
        It has <NUM> tags (maybe).
        """)
    subst_widget_four = dedent("""\
        This widget is a two-line string.
        It has four tags (maybe).
        """)

    dict_widget_four = {'NOUN' : 'widget',
                        'NUM' : 'four'}

    def test_Utils_TemplateSubst_GoodSubst(self):
        from opan.utils import template_subst as ts

        self.assertEqual(self.subst_widget_four,
                ts(self.template, self.dict_widget_four))


class TestOpanUtilsBaseAssertNPFArray(unittest.TestCase):

    from opan.error import OpanError
    from opan.test.utils import assertErrorAndTypecode

    class TestError(OpanError):
        NONE = 'NONE'
        NOT_FLOAT = 'NOT_FLOAT'
        NOT_ARRAY = 'NOT_ARRAY'
        NO_MEMBER = 'NO_MEMBER'

    def test_Utils_AssertNPFArray_TestObjIs1DFloatArray(self):
        import numpy as np
        from opan.utils import assert_npfloatarray as a_npfa

        try:
            a_npfa(np.float_(np.array(range(5))), None, "1D Float array",
                    self.TestError, self.TestError.NONE,
                    "No actual error; this should never be raised")
        except Exception:
            self.fail("Assertion failed on valid 1-D ndarray")

    def test_Utils_AssertNPFArray_TestObjIs2DFloatArray(self):
        import numpy as np
        from opan.utils import assert_npfloatarray as a_npfa

        try:
            a_npfa(np.float_(np.array(range(6)).reshape((2,3))),
                    None, "1D float array",
                    self.TestError, self.TestError.NONE,
                    "No actual error; this should never be raised")
        except Exception:
            self.fail("Assertion failed on valid 2-D ndarray")

    def test_Utils_AssertNPFArray_TestObjArrayNotFloat(self):
        import numpy as np
        from opan.utils import assert_npfloatarray as a_npfa

        self.assertErrorAndTypecode(self.TestError, a_npfa,
                self.TestError.NOT_FLOAT,
                np.array(range(5)), None, "1D int array",
                self.TestError, self.TestError.NOT_FLOAT,
                "ASSERTION FAIL: Non-float array")

    def test_Utils_AssertNPFArray_TestObjNotArray(self):
        import numpy as np
        from opan.utils import assert_npfloatarray as a_npfa

        # NumPy numeric type, but not an array
        self.assertErrorAndTypecode(self.TestError, a_npfa,
                self.TestError.NOT_ARRAY,
                np.float_(5.5), None, "Bare float",
                self.TestError, self.TestError.NOT_ARRAY,
                "ASSERTION FAIL: Non-array")

        # Non-ndarray type
        self.assertErrorAndTypecode(self.TestError, a_npfa,
                self.TestError.NOT_ARRAY,
                self.TestError, None, "Bare float",
                self.TestError, self.TestError.NOT_ARRAY,
                "ASSERTION FAIL: Non-array")

    def test_Utils_AssertNPFArray_TestVarIsFloatArray(self):
        import numpy as np
        from opan.utils import assert_npfloatarray as a_npfa

        class TestClass(object):
            import numpy as np
            testvar = np.float_(np.array(range(5)))

        try:
            a_npfa(TestClass, 'testvar', "Array member",
                    self.TestError, self.TestError.NONE,
                    "No actual error; this should never be raised")
        except Exception:
            self.fail("Assertion failed on a valid ndarray member")

    def test_Utils_AssertNPFArray_TestVarNotPresent(self):
        import numpy as np
        from opan.utils import assert_npfloatarray as a_npfa

        self.assertErrorAndTypecode(self.TestError, a_npfa,
                self.TestError.NO_MEMBER,
                self.TestError, 'not_there', "Absent member",
                self.TestError, self.TestError.NO_MEMBER,
                "ASSERTION FAIL: Member not found")



def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanUtilsBaseMisc),
                tl.loadTestsFromTestCase(TestOpanUtilsBaseCheckGeom),
                tl.loadTestsFromTestCase(TestOpanUtilsBaseTemplateSubst),
                tl.loadTestsFromTestCase(TestOpanUtilsBaseAssertNPFArray)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

