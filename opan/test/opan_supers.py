#-------------------------------------------------------------------------------
# Name:        opan_supers
# Purpose:     High-level test superclasses
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     12 Mar 2016
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


class TestSuperOpanGrad(unittest.TestCase):

    def test_SuperOpanGradIsAbstract(self):
        from opan.grad import SuperOpanGrad

        self.assertRaises(NotImplementedError, SuperOpanGrad)

## end class TestSuperOpanGrad


class TestSuperOpanHess(unittest.TestCase):

    def test_SuperOpanHessIsAbstract(self):
        from opan.hess import SuperOpanHess

        self.assertRaises(NotImplementedError, SuperOpanHess)

## end class TestSuperOpanHess



class SuperOrca(unittest.TestCase):
    # Superclass for all ORCA test case superclasses

    # Imports
    import os

    # Constants
    testdir = 'orca_test_dir'
    resourcedir = os.path.join('test','resource','orca')

## end class SuperOrca


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestSuperOpanGrad),
                tl.loadTestsFromTestCase(TestSuperOpanHess)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")


