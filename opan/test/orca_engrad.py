#-------------------------------------------------------------------------------
# Name:        orca_engrad
# Purpose:     Tests for opan.grad.OrcaEngrad
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     12 Mar 2016
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
from opan.test.opan_supers import SuperOrca


# ============================  OrcaEngrad ================================= #

class SuperOrcaEngrad(SuperOrca):
    # Superclass for all ORCA engrad test cases

    # Imports
    from textwrap import dedent
    import numpy as np
    from opan.test.utils import assertErrorAndTypecode

    # Superclass constants
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
    class names(object):
        numats = 'numats'
        energy = 'energy'
        grad = 'grad'
        geom = 'geom'
        atomicnum = 'atomicnum'
        E = frozenset([numats, energy, grad, geom, atomicnum])

    bad_block_substs = {
            names.numats: ('umber of', 'asdlkjf'),
            names.energy: ('urrent total en', 'alksdjfdsflkj'),
            names.grad: ('Eh/bohr', 'asdlkjsdf'),
            names.geom: ('inates in Bohr', 'asldkjfas;ldkfj')
                        }

    trunc_block_substs = {
            names.grad: ('-0.000007912155', '#-0.000007912155'),
            names.geom: ('1     2.6798241', '#1     2.6798241')
                            }

    bad_data_substs = {
            names.numats: ('of atoms\n#\n 4', 'of atoms\n#\n 8'),
            names.atomicnum: (' 29    -1.5545432', '229    -1.5545432')
                        }

    atoms = ['CU', 'O', 'H', 'H']
    gradient = np.array([0.000004637000, 0.000001922807, 0.000002366827,
                0.000010704036, -0.000004735732, -0.000005580943,
                0.000011997130, -0.000008391237, -0.000007912155,
                0.000011493781, -0.000008177479, -0.000008233658])
    energy = -1715.759691151236
    geom = np.array([-1.5545432, -0.4923021, -0.4922193,
                            1.8640436, 0.3660366, 0.3660223,
                            2.6798241, 1.9892213, -0.1622520,
                            2.6798521, -0.1622023, 1.9892043])

## end class SuperOrcaEngrad


class TestOrcaEngradKnownGood(SuperOrcaEngrad):
    # Testing errors related to file parsing of a known-good ENGRAD file

    @classmethod
    def setUpClass(cls):
        from opan.test.utils import setUpTestDir

        # Set up the directory
        setUpTestDir(cls.testdir)

        # Write the file
        with open(cls.file_name, 'w') as f:
            f.write(cls.file_text_good)


    @classmethod
    def tearDownClass(cls):
        import os
        from opan.test.utils import tearDownTestDir

        # Delete the engrad file
        os.remove(cls.file_name)

        # Remove the working directory
        tearDownTestDir(cls.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.grad import OrcaEngrad

        # Create the object
        self.oe = OrcaEngrad(path=self.file_name)

        # Enable long messages
        self.longMessage = True

    def test_ENGRAD_KnownGoodAtomVec(self):
        # Confirm the values coming out of the ENGRAD match the known-good
        #  example file.

        # Confirm the atom list is good
        for i in range(len(self.oe.atom_syms)):
            self.assertEqual(self.oe.atom_syms[i], self.atoms[i])

    def test_ENGRAD_KnownGoodGradient(self):
        # Confirm the known-good gradient matches what's expected.

        # Confirm the gradient vector is good
        for i in range(self.oe.gradient.shape[0]):
            self.assertAlmostEqual(self.oe.gradient[i],
                        self.gradient[i], delta=1e-10,
                        msg="Gradient index " + str(i))

    def test_ENGRAD_KnownGoodEnergy(self):
        # Confirm the known-good energy matches what's expected
        self.assertAlmostEqual(self.energy, self.oe.energy, delta=1e-12)

    def test_ENGRAD_KnownGoodGeom(self):
        for i in range(self.oe.geom.shape[0]):
            self.assertAlmostEqual(self.oe.geom[i],
                        self.geom[i], delta=1e-7,
                        msg="Coordinate index " + str(i))

    def test_ENGRAD_KnownGoodCheckGeomMatches(self):
        self.assertTrue(self.oe.check_geom(self.oe.geom, self.oe.atom_syms))

## end class TestOrcaEngradKnownGood


class TestOrcaEngradMissingBlocks(SuperOrcaEngrad):
    # Ensuring importing a non-ENGRAD file throws the right errors

    @classmethod
    def setUpClass(cls):
        # Set up the directory and add munged files
        from opan.test.utils import setUpTestDir

        setUpTestDir(cls.testdir)

        # Write the files
        for bname in cls.bad_block_substs.keys():
            with open(cls.file_name + bname, 'w') as f:
                f.write(cls.file_text_good
                                    .replace(*cls.bad_block_substs[bname]))

    @classmethod
    def tearDownClass(cls):
        # Remove any engrad files and try to remove the temp directory

        import os
        from opan.test.utils import tearDownTestDir

        # Try to remove the files
        [os.remove(cls.file_name + bname) for bname in
                                            cls.bad_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(cls.testdir)

    def test_ENGRAD_MissingBlockNumAts(self):

        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                    GradError.NUMATS, path=(self.file_name + self.names.numats))

    def test_ENGRAD_MissingBlockEnergy(self):

        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                    GradError.ENERGY, path=(self.file_name + self.names.energy))

    def test_ENGRAD_MissingBlockGrad(self):

        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                GradError.GRADBLOCK, path=(self.file_name + self.names.grad))

    def test_ENGRAD_MissingBlockGeom(self):

        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                GradError.GEOMBLOCK, path=(self.file_name + self.names.geom))

## end def TestOrcaEngradMissingBlocks


class TestOrcaEngradTruncatedBlocks(SuperOrcaEngrad):
    # Ensuring importing a file with incomplete grad or geom block throws
    #  the right errors

    @classmethod
    def setUpClass(cls):
        # Set up the directory and add munged files
        from opan.test.utils import setUpTestDir

        setUpTestDir(cls.testdir)

        # Write the files
        for bname in cls.trunc_block_substs.keys():
            with open(cls.file_name + bname, 'w') as f:
                f.write(cls.file_text_good
                                    .replace(*cls.trunc_block_substs[bname]))


    @classmethod
    def tearDownClass(cls):
        # Remove any engrad files and try to remove the temp directory

        import os
        from opan.test.utils import tearDownTestDir

        # Try to remove the files
        [os.remove(cls.file_name + bname) for bname in
                                        cls.trunc_block_substs.keys()]

        # Try to remove the directory
        tearDownTestDir(cls.testdir)

    def test_ENGRAD_TruncatedBlockGrad(self):

        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                GradError.GRADBLOCK, path=(self.file_name + self.names.grad))

    def test_ENGRAD_TruncatedBlockGeom(self):

        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                GradError.GEOMBLOCK, path=(self.file_name + self.names.geom))

## end class TestOrcaEngradTruncatedBlocks


class TestOrcaEngradBadData(SuperOrcaEngrad):
    # Ensuring files with valid formatting but invalid data raise errors

    @classmethod
    def setUpClass(cls):
        # Set up the directory and add munged files
        from opan.test.utils import setUpTestDir

        # Create the test directory
        setUpTestDir(cls.testdir)

        # Write the files
        for dname in cls.bad_data_substs.keys():
            with open(cls.file_name + dname, 'w') as f:
                f.write(cls.file_text_good
                                    .replace(*cls.bad_data_substs[dname]))

    @classmethod
    def tearDownClass(cls):
        # Tear down test directory and remove files

        import os
        from opan.test.utils import tearDownTestDir

        # Try to remove the files
        [os.remove(cls.file_name + dname) for dname in
                                            cls.bad_data_substs.keys()]

        # Remove the test directory
        tearDownTestDir(cls.testdir)


    def test_ENGRAD_BadDataAtomicNum(self):
        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                    GradError.GEOMBLOCK, path=(self.file_name
                                        + self.names.atomicnum))

    def test_ENGRAD_BadDataNumAtoms(self):
        from opan.error import GradError
        from opan.grad import OrcaEngrad

        self.assertErrorAndTypecode(GradError, OrcaEngrad,
                    GradError.GRADBLOCK, path=(self.file_name
                                        + self.names.numats))

## end class TestOrcaEngradBadData


class TestOrcaEngradBadUsage(SuperOrcaEngrad):

    @classmethod
    def setUpClass(cls):
        from opan.test.utils import setUpTestDir

        # Set up the directory
        setUpTestDir(cls.testdir)

        # Write the file
        with open(cls.file_name, 'w') as f:
            f.write(cls.file_text_good)


    @classmethod
    def tearDownClass(cls):
        import os
        from opan.test.utils import tearDownTestDir

        # Delete the engrad file
        os.remove(cls.file_name)

        # Remove the working directory
        tearDownTestDir(cls.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.grad import OrcaEngrad

        # Create the object
        self.oe = OrcaEngrad(path=self.file_name)

        # Enable long messages
        self.longMessage = True

    def test_ENGRAD_BadUsageOverwrite(self):
        from opan.error import GradError
        self.assertErrorAndTypecode(GradError, self.oe.load,
                                GradError.OVERWRITE, path=self.oe.engrad_path)

## end class TestOrcaEngradBadUsage


class TestOrcaEngradLiveData(SuperOrcaEngrad):

    def test_ENGRAD_LiveData(self):

        import os
        from opan.grad import OrcaEngrad
        from opan.error import GradError

        for fname in os.listdir(self.resourcedir):
            if fname.startswith("test_orca") and fname.endswith("engrad"):
                print("\nTesting file '" + fname + "' ... ")
                try:
                    OrcaEngrad(path=os.path.join(self.resourcedir, fname))
                except (IOError, GradError) as e: # pragma: no cover
                    self.longMessage = True
                    self.fail("Load of test file '" + str(fname) +
                            "' failed:\n" + str(e))

## end class TestOrcaEngradLiveData


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOrcaEngradBadData),
                tl.loadTestsFromTestCase(TestOrcaEngradBadUsage),
                tl.loadTestsFromTestCase(TestOrcaEngradKnownGood),
                tl.loadTestsFromTestCase(TestOrcaEngradLiveData),
                tl.loadTestsFromTestCase(TestOrcaEngradMissingBlocks),
                tl.loadTestsFromTestCase(TestOrcaEngradTruncatedBlocks)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")
