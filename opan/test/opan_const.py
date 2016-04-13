#-------------------------------------------------------------------------------
# Name:        opan_const
# Purpose:     Test objects for opan.const
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     10 Mar 2016
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


class TestOpanEnumValueCheck(unittest.TestCase):

    def test_OpanEnum_ValueCheck(self):
        from opan.const import EnumDispDirection as EDD

        # Representative value in a representative Enum
        self.assertTrue(EDD.NEGATIVE in EDD)

    def test_OpanEnum_IterCheck(self):
        from opan.const import EnumDispDirection as EDD

        self.assertSetEqual({'NEGATIVE', 'NO_DISP', 'POSITIVE'},
                            set(k for k in EDD))


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanEnumValueCheck)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

