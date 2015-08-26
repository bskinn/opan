#-------------------------------------------------------------------------------
# Name:        tests
# Purpose:     Definitions of all unit tests for the OpenAnharmonic package.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     30 Jul 2015
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

# Module-level imports
import unittest


# ============================  ORCA_ENGRAD ================================= #

class SuperORCAEngrad(unittest.TestCase):
    # Superclass for all engrad test cases

    # Imports
    from textwrap import dedent
    import numpy as np

    # Class variables
    testdir = 'orca_test_dir'
    file_name = 'test.engrad'
    file_text_good = dedent("""\
        #
        # Number of atoms
        #
         4
        #
        # The current total energy in Eh
        #
          -1715.759691151236
        #
        # The current gradient in Eh/bohr
        #
               0.000004637000
               0.000001922807
               0.000002366827
               0.000010704036
              -0.000004735732
              -0.000005580943
               0.000011997130
              -0.000008391237
              -0.000007912155
               0.000011493781
              -0.000008177479
              -0.000008233658
        #
        # The atomic numbers and current coordinates in Bohr
        #
          29    -1.5545432   -0.4923021   -0.4922193
           8     1.8640436    0.3660366    0.3660223
           1     2.6798241    1.9892213   -0.1622520
           1     2.6798521   -0.1622023    1.9892043
        """)
    class blocknames(object):
        numats = 'numats'
        energy = 'energy'
        grad = 'grad'
        geom = 'geom'
        E = frozenset([numats, energy, grad, geom])

    bad_block_substs = {
            'numats': ('umber of', 'asdlkjf'),
            'energy': ('urrent total en', 'alksdjfdsflkj'),
            'grad': ('urrent gradient', 'asdlkjsdf'),
            'geom': ('atomic numbers and current', 'asldkjfas;ldkfj')
                        }

    atoms = np.array([['CU'],['O'],['H'],['H']])
    gradient = np.matrix([0.000004637000, 0.000001922807, 0.000002366827, \
                0.000010704036, -0.000004735732, -0.000005580943, \
                0.000011997130, -0.000008391237, -0.000007912155, \
                0.000011493781, -0.000008177479, -0.000008233658]).transpose()
    energy = -1715.759691151236
    coords = np.matrix([-1.5545432, -0.4923021, -0.4922193, \
                            1.8640436, 0.3660366, 0.3660223, \
                            2.6798241, 1.9892213, -0.1622520, \
                            2.6798521, -0.1622023, 1.9892043]).transpose()


class TestORCAEngradKnownGood(SuperORCAEngrad):
    # Testing errors related to file parsing of a known-good ENGRAD file

    @classmethod
    def setUpClass(self):
        import os

        # Check if test directory already exists (or file of same name);
        #  error if so
        if os.path.isdir(self.testdir) or os.path.isfile(self.testdir):
            raise(IOError("Cannot create new test directory!"))

        # Create and change to test directory
        os.mkdir(self.testdir)
        os.chdir(self.testdir)

        # Write the file
        with open(self.file_name, 'w') as f:
            f.write(self.file_text_good)


    @classmethod
    def tearDownClass(self):
        import os

        # Delete the engrad file
        os.remove(self.file_name)

        # Switch to parent directory
        os.chdir(os.path.pardir)

        # Try to remove the temp directory
        os.rmdir(self.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.grad import ORCA_ENGRAD as OE

        # Create the object
        self.oe = OE(self.file_name)

    def test_ENGRAD_KnownGoodAtomVec(self):
        # Confirm the values coming out of the ENGRAD match the known-good
        #  example file included as the below string.

        # Confirm the atom list is good
        for i in range(4):
            self.assertEqual(self.oe.atom_syms[i,0], self.atoms[i,0])

    def test_ENGRAD_KnownGoodGradient(self):
        # Confirm the known-good gradient matches what's expected.

        # Confirm the gradient vector is good
        self.longMessage = True
        for i in range(12):
            self.assertAlmostEqual(self.oe.gradient[i,0], \
                        self.gradient[i,0], delta=1e-10, \
                        msg="Gradient index " + str(i))

    def test_ENGRAD_KnownGoodEnergy(self):
        # Confirm the known-good energy matches what's expected
        self.assertAlmostEqual(self.energy, self.oe.energy, delta=1e-12)

    def test_ENGRAD_KnownGoodGeom(self):
        self.longMessage = True
        for i in range(12):
            self.assertAlmostEqual(self.oe.geom_vec[i,0], \
                        self.coords[i,0], delta=1e-7, \
                        msg="Coordinate index " + str(i))

    def test_ENGRAD_KnownGoodInitFlagDefined(self):
        self.assertTrue('initialized' in self.oe.__dict__)
        self.assertTrue(self.oe.initialized)


## end class TestORCAEngradKnownGood


class TestORCAEngradMissingBlocks(SuperORCAEngrad):
    # Ensuring importing a non-ENGRAD file throws the right errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        import os

        # Check if test directory already exists (or file of same name);
        #  error if so
        if os.path.isdir(self.testdir) or os.path.isfile(self.testdir):
            raise(IOError("Cannot create new test directory!"))

        # Create and change to test directory
        os.mkdir(self.testdir)
        os.chdir(self.testdir)

        # Write the files
        for bname in self.blocknames.E:
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good.replace \
                                            (*self.bad_block_substs[bname]))


    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in self.blocknames.E]

        # Switch to parent directory
        os.chdir(os.path.pardir)

        # Try to remove the temp directory
        os.rmdir(self.testdir)

    def test_ENGRAD_MissingBlockNumAts(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        # Munged number-of-atoms block
        self.assertRaises(GRADError, ORCA_ENGRAD, \
                                    self.file_name + self.blocknames.numats)
        # Ensure correct typecode; suppress repeated error
        try:
            oe = ORCA_ENGRAD(self.file_name + self.blocknames.numats)
        except GRADError as ge:
            self.assertEqual(ge.tc, GRADError.numats)
        except:
            pass


## end def TestORCAEngradMungeGrad


# ============================  ORCA_ERROR ================================== #

class TestOPANErrorInitErrors(unittest.TestCase):
    # Testing errors that should be thrown on initialization

    def test_ORCAError_init_NotImplemented(self):
        # Must import
        from opan.error import OPANError

        # Confirm OPANError parent class as abstract
        self.assertRaises(NotImplementedError, OPANError, \
                    "tc", "msg", "src")

    def test_XYZError_init_BadTypecode(self):
        # Must import
        from opan.error import XYZError

        # Confirm KeyError raised when invalid typecode passed
        self.assertRaises(KeyError, XYZError, \
                    "INVALID TYPECODE", "msg", "src")

class TestXYZErrorInitConfig(unittest.TestCase):
    # XYZError used as a representative OPANError subclass

    # Class-level constants
    tc = 'nonprl'
    msg = "Test message"
    src = "Test source"
    subclass_name = "XYZError"

    def test_XYZError_init_SubclassName(self):
        # Confirm subclass name is retrieved correctly
        from opan.error import XYZError as XE
        self.assertEqual(XE(self.tc, self.msg, self.src).subclass_name, \
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



if __name__ == '__main__':
    unittest.main(verbosity=2)
