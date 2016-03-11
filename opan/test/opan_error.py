#-------------------------------------------------------------------------------
# Name:        opan_error
# Purpose:     Test objects for opan.error
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     10 Mar 2016
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


class TestOpanErrorTypecodeCheck(unittest.TestCase):

    def test_OpanError_TypecodeCheck(self):
        from opan.error import XYZError

        # Representative typecode in a representative error class
        self.assertTrue(XYZError.XYZFILE in XYZError)


class TestOpanErrorInitErrors(unittest.TestCase):
    # Testing errors that should be thrown on initialization

    def test_OpanError_init_NotImplemented(self):
        # Must import
        from opan.error import OpanError

        # Confirm OpanError parent class as abstract
        self.assertRaises(NotImplementedError, OpanError, "tc", "msg", "src")

    def test_XYZError_init_BadTypecode(self):
        # Must import
        from opan.error import XYZError

        # Confirm KeyError raised when invalid typecode passed
        self.assertRaises(KeyError, XYZError, "INVALID TYPECODE", "msg", "src")


class TestXYZErrorInitConfig(unittest.TestCase):
    # XYZError used as a representative OpanError subclass

    # Class-level constants
    tc = 'NONPRL'
    msg = "Test message"
    src = "Test source"
    subclass_name = "XYZError"

    def test_XYZError_init_SubclassName(self):
        # Confirm subclass name is retrieved correctly
        from opan.error import XYZError as XE
        self.assertEqual(XE(self.tc, self.msg, self.src).subclass_name,
                                                self.subclass_name)

    def test_XYZError_init_TypecodeStored(self):
        # Confirm typecode name stored correctly
        from opan.error import XYZError as XE
        self.assertEqual(XE(self.tc, self.msg, self.src).tc, self.tc)

    def test_XYZError_init_MessageStrStored(self):
        # Confirm message string stored correctly
        from opan.error import XYZError as XE
        self.assertEqual(XE(self.tc, self.msg, self.src).msg, self.msg)

    def test_XYZError_init_SourceStrStored(self):
        # Confirm source string stored correctly
        from opan.error import XYZError as XE
        self.assertEqual(XE(self.tc, self.msg, self.src).src, self.src)


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanErrorInitErrors),
                tl.loadTestsFromTestCase(TestOpanErrorTypecodeCheck),
                tl.loadTestsFromTestCase(TestXYZErrorInitConfig)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")


