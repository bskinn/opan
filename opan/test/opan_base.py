#-------------------------------------------------------------------------------
# Name:        opan_base
# Purpose:     Base-level tests for opan
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     20 Jun 2016
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


class TestOpanBaseImports(unittest.TestCase):

    def test_OpanBaseImportNumPy(self):
        try:
            import numpy
        except Exception:
            self.fail(msg="NumPy import failed.")

    def test_OpanBaseImportSciPy(self):
        try:
            from scipy import linalg
            from scipy import sparse
        except Exception:
            self.fail(msg="SciPy subpackage import(s) failed.")

    def test_OpanBaseImportH5Py(self):
        try:
            import h5py
        except Exception:
            self.fail(msg="H5Py import failed.")

    def test_OpanBaseImportOpan(self):
        try:
            import opan
        except Exception:
            self.fail(msg="opan import failed.")


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanBaseImports),
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")


