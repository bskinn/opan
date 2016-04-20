#-------------------------------------------------------------------------------
# Name:        opan_utils_vector
# Purpose:     Test objects for opan.utils.vector
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     19 Apr 2016
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


class TestOpanUtilsVectorMisc(unittest.TestCase):

    def test_Utils_Vector_ParallelCheck_Good(self):
        from opan.utils.vector import parallel_check as pc
        import numpy as np

        # Parallel vectors
        self.assertTrue(pc(np.array([1, 2, 3]),
                           np.array([1.5, 3, 4.5])))

        # Anti-parallel vectors
        self.assertTrue(pc(np.array([-1, 5.3, 3.1]),
                           np.array([3, -15.9, -9.3])))

        # Non-(anti-)parallel vectors
        self.assertFalse(pc(np.array([4.8, 0.35, -1.822]),
                            np.array([-1.3, 3.77, 19.14])))

    def test_Utils_Vector_ParallelCheck_BadShape(self):
        from opan.utils.vector import parallel_check as pc
        import numpy as np

        self.assertRaises(ValueError, pc,
                          np.array([[1, 2, 3], [3, 2, 1]]),
                          np.array([2, 4, 9]))

    def test_Utils_Vector_ParallelCheck_LenMismatch(self):
        from opan.utils.vector import parallel_check as pc
        import numpy as np

        self.assertRaises(ValueError, pc,
                          np.array(range(3)),
                          np.array(range(4)))


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanUtilsVectorMisc)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")
