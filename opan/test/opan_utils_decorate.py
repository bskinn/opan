#-------------------------------------------------------------------------------
# Name:        opan_utils_decorate
# Purpose:     Test objects for opan.utils.decorate
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     12 Apr 2016
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


class TestOpanUtilsDecorate(unittest.TestCase):

    from opan.utils.decorate import arraysqueeze as arsq

    @classmethod
    def setUpClass(cls):
        cls.longMessage = True

    @staticmethod
    @arsq(0)
    def fxn_2p_1d(p1, p2):
        return (p1, p2)

    @staticmethod
    @arsq(0,1)
    def fxn_2p_2d(p1, p2):
        return (p1, p2)

    @staticmethod
    @arsq(1, 'test')
    def fxn_2p_1d_k(p1, p2, **kwargs):
        return p1, p2, (kwargs['test'] if 'test' in kwargs else None)


    def test_arraysqueeze_GoodPartialWrap(self):
        import numpy as np

        # Store the result for detailed checking
        ret = self.fxn_2p_1d(1,2)

        # First one needs to be an ndarray
        self.assertIsInstance(ret[0], np.ndarray, 
                msg="First argument not converted to ndarray")
        self.assertTupleEqual(ret[0].shape, (1,),
                msg="Bad ndarray size for wrapped first argument")
        self.assertEqual(ret[0][0], 1,
                msg="Corrupted value in wrapped first argument")

        # Second one needs to still be an int
        self.assertIsInstance(ret[1], int,
                msg="Second argument not left as int")
        self.assertEqual(ret[1], 2,
                msg="Corrupted value in unmodified second argument")

    def test_arraysqueeze_GoodTwoArgWrap(self):
        import numpy as np

        # Store the result for detailed checking
        ret = self.fxn_2p_2d(1,2)

        # Both need to be ndarrays
        self.assertIsInstance(ret[0], np.ndarray,
                msg="First argument not converted to ndarray")
        self.assertTupleEqual(ret[0].shape, (1,),
                msg="Bad ndarray size for wrapped first argument")
        self.assertEqual(ret[0][0], 1,
                msg="Corrupted value in wrapped first argument")

        self.assertIsInstance(ret[1], np.ndarray,
                msg="Second argument not converted to ndarray")
        self.assertTupleEqual(ret[1].shape, (1,),
                msg="Bad ndarray size for wrapped second argument")
        self.assertEqual(ret[1][0], 2,
                msg="Corrupted value in wrapped second argument")

    def test_arraysqueeze_GoodArrayArgWrap(self):
        import numpy as np
        from opan.utils import pack_tups

        # Store the array result for detailed checking
        ary = np.array(range(1,9,2)).reshape((4,1))
        ret = self.fxn_2p_1d(ary,"x")

        # Confirm converted array matches expectations
        for i,t in enumerate(pack_tups(ary, ret[0])):
            self.assertEqual(t[0][0], t[1],
                msg="Mismatch at position {0} (Pre: {1}; Post: {2})"
                    .format(i, t[0][0], t[1]))

    def test_arraysqueeze_GoodKeywordArgWrap(self):
        import numpy as np

        # Result function with keyword
        ret = self.fxn_2p_1d_k(1, 2, test=3)

        # Confirm that the parameters are wrapped properly
        self.assertIsInstance(ret[1], np.ndarray,
                msg="Second argument not converted to ndarray")
        self.assertIsNotNone(ret[2],
                msg="Keyword argument improperly returned as None")
        self.assertIsInstance(ret[2], np.ndarray,
                msg="Keyword argument not converted to ndarray")

    def test_arraysqueeze_BadInitParam(self):

        with self.assertRaises(ValueError,
                    msg="Error not raised on float param index"):
            @self.arsq(1.5)
            def test_fxn(p1, p2):
                pass

        with self.assertRaises(ValueError,
                    msg="Error not raised on object param index"):
            @self.arsq(ValueError)
            def test_fxn(p1, p2):
                pass

## end class TestOpanUtilsDecorate


def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(TestOpanUtilsDecorate)
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

