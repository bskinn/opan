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


class TestOpanUtilsBase(unittest.TestCase):

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


def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(TestOpanUtilsBase)
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

