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


# =============================   ORCA   ==================================== #

class SuperORCA(unittest.TestCase):
    # Superclass for all ORCA test case superclasses

    # Imports
    import os

    # Constants
    testdir = 'orca_test_dir'
    resourcedir = os.path.join('resource','orca')

## end class SuperORCA


# ============================  ORCA_ENGRAD ================================= #

class SuperORCAEngrad(SuperORCA):
    # Superclass for all ORCA engrad test cases

    # Imports
    from textwrap import dedent
    import numpy as np

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

## end class SuperORCAEngrad


class TestORCAEngradKnownGood(SuperORCAEngrad):
    # Testing errors related to file parsing of a known-good ENGRAD file

    @classmethod
    def setUpClass(self):

        # Set up the directory
        setUpTestDir(self.testdir)

        # Write the file
        with open(self.file_name, 'w') as f:
            f.write(self.file_text_good)


    @classmethod
    def tearDownClass(self):
        import os

        # Delete the engrad file
        os.remove(self.file_name)

        # Remove the working directory
        tearDownTestDir(self.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.grad import ORCA_ENGRAD

        # Create the object
        self.oe = ORCA_ENGRAD(self.file_name)

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
        for i in range(self.oe.geom_vec.shape[0]):
            self.assertAlmostEqual(self.oe.geom_vec[i],
                        self.geom[i], delta=1e-7,
                        msg="Coordinate index " + str(i))

    def test_ENGRAD_KnownGoodInitFlagDefined(self):
        self.assertTrue('initialized' in self.oe.__dict__)
        if 'initialized' in self.oe.__dict__:
            self.assertTrue(self.oe.initialized)

## end class TestORCAEngradKnownGood


class TestORCAEngradMissingBlocks(SuperORCAEngrad):
    # Ensuring importing a non-ENGRAD file throws the right errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.bad_block_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.bad_block_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in
                                            self.bad_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_ENGRAD_MissingBlockNumAts(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.numats, self.file_name + self.names.numats)

    def test_ENGRAD_MissingBlockEnergy(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.en, self.file_name + self.names.energy)

    def test_ENGRAD_MissingBlockGrad(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.gradblock, self.file_name + self.names.grad)

    def test_ENGRAD_MissingBlockGeom(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.geomblock, self.file_name + self.names.geom)

## end def TestORCAEngradMissingBlocks


class TestORCAEngradTruncatedBlocks(SuperORCAEngrad):
    # Ensuring importing a file with incomplete grad or geom block throws
    #  the right errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.trunc_block_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.trunc_block_substs[bname]))


    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in
                                        self.trunc_block_substs.keys()]

        # Try to remove the directory
        tearDownTestDir(self.testdir)

    def test_ENGRAD_TruncatedBlockGrad(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.gradblock, self.file_name + self.names.grad)

    def test_ENGRAD_TruncatedBlockGeom(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.geomblock, self.file_name + self.names.geom)

## end class TestORCAEngradTruncatedBlocks


class TestORCAEngradBadData(SuperORCAEngrad):
    # Ensuring files with valid formatting but invalid data raise errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        # Create the test directory
        setUpTestDir(self.testdir)

        # Write the files
        for dname in self.bad_data_substs.keys():
            with open(self.file_name + dname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.bad_data_substs[dname]))

    @classmethod
    def tearDownClass(self):
        # Tear down test directory and remove files

        import os

        # Try to remove the files
        [os.remove(self.file_name + dname) for dname in
                                            self.bad_data_substs.keys()]

        # Remove the test directory
        tearDownTestDir(self.testdir)


    def test_ENGRAD_BadDataAtomicNum(self):
        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.geomblock, self.file_name
                                        + self.names.atomicnum)

    def test_ENGRAD_BadDataNumAtoms(self):
        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD,
                    GRADError.gradblock, self.file_name
                                        + self.names.numats)

## end class TestORCAEngradBadData


class TestORCAEngradLiveData(SuperORCAEngrad):

    def test_ENGRAD_LiveData(self):

        import os
        from opan.grad import ORCA_ENGRAD
        from opan.error import GRADError

        for fname in os.listdir(self.resourcedir):
            if fname[:9] == "test_orca" and fname[-6:] == "engrad":
                print("\nTesting file '" + fname + "' ... ")
                try:
                    ORCA_ENGRAD(os.path.join(self.resourcedir, fname))
                except (IOError, GRADError) as e: # pragma: no cover
                    self.longMessage = True
                    self.fail("Load of test file '" + str(fname) +
                            "' failed:\n" + str(e))

## end class TestORCAEngradLiveData


# ============================  OPAN_ERROR ================================== #

class TestOPANErrorInitErrors(unittest.TestCase):
    # Testing errors that should be thrown on initialization

    def test_OPANError_init_NotImplemented(self):
        # Must import
        from opan.error import OPANError

        # Confirm OPANError parent class as abstract
        self.assertRaises(NotImplementedError, OPANError, "tc", "msg", "src")

    def test_XYZError_init_BadTypecode(self):
        # Must import
        from opan.error import XYZError

        # Confirm KeyError raised when invalid typecode passed
        self.assertRaises(KeyError, XYZError, "INVALID TYPECODE", "msg", "src")

## end class TestOPANErrorInitErrors


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

## end class TestXYZErrorInitConfig


# =============================  ORCA_HESS  ================================= #

class SuperORCAHess(SuperORCA):
    # Superclass for all ORCA .hess test cases

    # Imports
    from textwrap import dedent
    import numpy as np

    # Superclass constants

    file_text_good = dedent("""\

    $orca_hessian_file

    $act_atom
      4

    $act_coord
      2

    $act_energy
          -40.451555

    $hessian
    15
                      0          1          2          3          4          5
          0       0.569133  -0.000007  -0.000000  -0.341891   0.000009   0.000000
          1      -0.000004   0.569121  -0.000000   0.000003  -0.042495  -0.000000
          2      -0.000000   0.000000   0.569119   0.000000   0.000000  -0.042494
          3      -0.341984  -0.000001   0.000000   0.364262  -0.000006   0.000000
          4       0.000004  -0.042439  -0.000000  -0.000002   0.041132   0.000000
          5       0.000000  -0.000000  -0.042433  -0.000000  -0.000000   0.041127
          6      -0.075716   0.094141  -0.000000  -0.007458   0.030752   0.000000
          7       0.094130  -0.308688   0.000000  -0.001778   0.002751  -0.000000
          8      -0.000000   0.000000  -0.042437  -0.000000  -0.000000  -0.001841
          9      -0.075717  -0.047067  -0.081528  -0.007456  -0.015377  -0.026632
         10      -0.047065  -0.108997  -0.115292   0.000889  -0.000694   0.001987
         11      -0.081519  -0.115289  -0.242125   0.001539   0.001987   0.001604
         12      -0.075717  -0.047067   0.081529  -0.007456  -0.015378   0.026632
         13      -0.047065  -0.108997   0.115292   0.000889  -0.000694  -0.001987
         14       0.081519   0.115290  -0.242125  -0.001539  -0.001987   0.001604
                      6          7          8          9         10         11
          0      -0.075684   0.094113   0.000000  -0.075685  -0.046981  -0.081374
          1       0.093944  -0.308614  -0.000000  -0.046971  -0.108902  -0.115083
          2      -0.000000  -0.000000  -0.042491  -0.081353  -0.115066  -0.241795
          3      -0.007493  -0.001776  -0.000000  -0.007493   0.000813   0.001407
          4       0.030795   0.002808  -0.000000  -0.015398  -0.000649   0.002005
          5       0.000000   0.000000  -0.001842  -0.026672   0.002005   0.001666
          6       0.076822  -0.101620  -0.000000   0.003178  -0.011625   0.001420
          7      -0.101487   0.328330  -0.000000   0.004579  -0.011317   0.001997
          8       0.000000   0.000000   0.041128   0.010779  -0.024479   0.001667
          9       0.003177   0.004642   0.010750   0.076820   0.050748   0.087909
         10      -0.011626  -0.011262  -0.024447   0.050743   0.112710   0.124317
         11       0.001422   0.002117   0.001603   0.087891   0.124307   0.256265
         12       0.003178   0.004642  -0.010750   0.003179   0.007046  -0.009361
         13      -0.011626  -0.011262   0.024447   0.007047   0.008159  -0.013236
         14      -0.001422  -0.002117   0.001602   0.009354   0.013234  -0.017804
                      12         13         14
          0      -0.075685  -0.046981   0.081374
          1      -0.046971  -0.108902   0.115083
          2       0.081353   0.115066  -0.241795
          3      -0.007493   0.000813  -0.001407
          4      -0.015398  -0.000649  -0.002005
          5       0.026672  -0.002005   0.001666
          6       0.003178  -0.011625  -0.001420
          7       0.004579  -0.011318  -0.001997
          8      -0.010779   0.024479   0.001667
          9       0.003179   0.007046   0.009361
         10       0.007047   0.008159   0.013236
         11      -0.009354  -0.013234  -0.017804
         12       0.076820   0.050748  -0.087909
         13       0.050743   0.112710  -0.124316
         14      -0.087891  -0.124307   0.256265

    $vibrational_frequencies
    15
        0        0.000000
        1        0.000000
        2        0.000000
        3        0.000000
        4        0.000000
        5        0.000000
        6     1292.217629
        7     1292.583415
        8     1292.939274
        9     1524.411513
       10     1524.709386
       11     3100.659679
       12     3239.096249
       13     3240.029637
       14     3241.018763

    $normal_modes
    15 15
                      0          1          2          3          4          5
          0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          2       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          3       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          4       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          5       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          6       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          7       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          8       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          9       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
         10       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
         11       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
         12       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
         13       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
         14       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
                      6          7          8          9         10         11
          0      -0.123657  -0.010576   0.000071  -0.000003   0.000001  -0.000165
          1      -0.010563   0.123653   0.000144  -0.000153  -0.000000  -0.000186
          2      -0.000083   0.000138  -0.124103  -0.000000  -0.000006   0.000000
          3      -0.097814  -0.008369   0.000056   0.000052   0.000001  -0.498497
          4       0.051382  -0.600708  -0.000702   0.500780   0.000180  -0.000010
          5       0.000401  -0.000672   0.603679   0.000179  -0.499975   0.000000
          6       0.542535  -0.174270  -0.000557   0.471703   0.000166   0.166288
          7       0.217937   0.039197  -0.000102   0.166582   0.000058  -0.470357
          8       0.000407  -0.000669   0.603630  -0.000178   0.500025   0.000000
          9       0.514246   0.154544  -0.191222  -0.236006   0.408156   0.167087
         10      -0.071907  -0.455650  -0.270601  -0.332665  -0.288800   0.236293
         11      -0.213339   0.253181   0.136154   0.288303   0.000116   0.409333
         12       0.514496   0.154117   0.190877  -0.235715  -0.408330   0.167087
         13      -0.071542  -0.456249   0.269684  -0.332870   0.288565   0.236294
         14       0.213521  -0.253482   0.135307  -0.288303  -0.000091  -0.409334
                      12         13         14
          0       0.000023   0.045812   0.080955
          1       0.000076   0.080969  -0.045802
          2      -0.093026   0.000078  -0.000018
          3      -0.000213  -0.425825  -0.750570
          4       0.000013   0.013314  -0.007506
          5      -0.015315   0.000013  -0.000003
          6       0.000207   0.199937  -0.207069
          7      -0.000562  -0.531035   0.616152
          8      -0.015328   0.000013  -0.000003
          9       0.238637  -0.160199  -0.003453
         10       0.337503  -0.223822  -0.031379
         11       0.569223  -0.410798  -0.041209
         12      -0.238904  -0.159800  -0.003544
         13      -0.337865  -0.223257  -0.031508
         14       0.569893   0.409845   0.041427

    #
    # The atoms: label  mass x y z
    #
    $atoms
    5
     C     12.0110     -0.000000     0.000000     0.000000
     H      1.0080      2.059801     0.000000     0.000000
     H      1.0080     -0.686600     1.942000     0.000000
     H      1.0080     -0.686600    -0.971000    -1.681821
     H      1.0080     -0.686600    -0.971000     1.681821

    $actual_temperature
      0.000000

    $dipole_derivatives
    15
        -0.030119    -0.000025     0.000000
         0.000010    -0.030119    -0.000000
        -0.000000     0.000000    -0.030159
        -0.118178     0.000019     0.000000
        -0.000029     0.072166     0.000000
        -0.000000     0.000000     0.072174
         0.050074     0.055974    -0.000000
         0.059226    -0.097241     0.000000
         0.000000    -0.000000     0.072183
         0.050094    -0.027985    -0.048461
        -0.028084     0.030186    -0.068789
        -0.048754    -0.068943    -0.049492
         0.050094    -0.027985     0.048461
        -0.028084     0.030186     0.068789
         0.048753     0.068944    -0.049492

    #
    # The IR spectrum
    #  wavenumber T**2 TX TY  TY
    #
    $ir_spectrum
    15
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
       1292.22      14.4422       3.7904       0.2737       0.0026
       1292.58      14.4744       0.3249      -3.7906      -0.0043
       1292.94      14.6749      -0.0022      -0.0044       3.8308
       1524.41       0.0002       0.0129      -0.0091       0.0000
       1524.71       0.0000      -0.0000       0.0000       0.0048
       3100.66       0.0203       0.0915       0.1093       0.0000
       3239.10      13.8923       0.0009       0.0032      -3.7272
       3240.03      14.9837       1.9037       3.3704       0.0031
       3241.02      16.1995       3.5028      -1.9825      -0.0007

    $polarizability_derivatives
    15
        -5.943413     2.971881     2.971618    -0.000247     0.000000     0.000000
        -0.000168    -4.203036     4.202646     2.965679     0.000000    -0.000000
        -0.000003    -0.000005     0.000002     0.000002     2.965322     4.193466
         9.078498     0.217638     0.217719    -0.000051     0.000001     0.000000
         0.000223     0.646863    -0.646912     1.996755    -0.000000     0.000001
         0.000006    -0.000007     0.000001     0.000003     1.996681    -0.640485
        -1.043996    -1.441613    -0.685449     2.585964     0.000000    -0.000001
         0.905935     8.074887    -0.011341    -2.036142    -0.000001    -0.000002
        -0.000023     0.000001    -0.000003    -0.000002    -1.269778     1.668773
        -1.044512    -0.874492    -1.252465    -1.292825    -2.239366    -0.330511
        -0.453132    -2.259789    -1.771464    -1.462185    -0.330621    -2.472069
        -0.784977    -1.016976    -5.965579    -0.330327    -1.844009    -2.608533
        -1.044517    -0.874493    -1.252466    -1.292826     2.239359     0.330504
        -0.453141    -2.259789    -1.771459    -1.462179     0.330622     2.472076
         0.784977     1.016966     5.965576     0.330327    -1.844010    -2.608531

    #
    # The Raman spectrum
    #  wavenumber Activity Depolarization
    #
    $raman_spectrum
    15
          0.00       0.0000       0.0000
          0.00       0.0000       0.0000
          0.00       0.0000       0.0000
          0.00       0.0000       0.0000
          0.00       0.0000       0.0000
          0.00       0.0000       0.0000
       1292.22       1.5128       0.7500
       1292.58       1.5204       0.7500
       1292.94       1.5313       0.7500
       1524.41      26.2882       0.7500
       1524.71      26.2096       0.7500
       3100.66     140.8157       0.0000
       3239.10      60.6507       0.7500
       3240.03      60.7749       0.7500
       3241.02      60.8277       0.7500

    $job_list
    15
        0    1    1    1
        1    1    1    1
        2    1    1    1
        3    1    1    1
        4    1    1    1

    $eigenvalues_mass_weighted_hessian
    15
        0          0.00000000000
        1          0.00000000000
        2          0.00000000000
        3          0.00000000000
        4          0.00000000000
        5          0.00000000000
        6          0.06319195757
        7          0.06322773793
        8          0.06326255701
        9          0.08794171404
       10          0.08797608543
       11          0.36383015346
       12          0.39704359307
       13          0.39727245265
       14          0.39751505083

    $eigenvectors_mass_weighted_hessian
    15 15
                           0                    1                    2                    3                    4                    5
          0            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          1            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          2            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          3            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          4            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          5            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          6            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          7            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          8            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
          9            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
         10            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
         11            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
         12            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
         13            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
         14            0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000        0.00000000000
                           6                    7                    8                    9                   10                   11
          0           -0.39494249803       -0.03377835105        0.00022673315       -0.00001016157        0.00000238076       -0.00056934730
          1           -0.03373814675        0.39492999897        0.00046131523       -0.00052911035       -0.00000098191       -0.00064301369
          2           -0.00026522927        0.00044031834       -0.39636743307       -0.00000021637       -0.00002182899        0.00000004376
          3           -0.09050154121       -0.00774361372        0.00005197726        0.00005230188        0.00000057151       -0.49849653717
          4            0.04754097370       -0.55580088222       -0.00064963081        0.50077959071        0.00017965318       -0.00000957941
          5            0.00037077591       -0.00062183177        0.55855083277        0.00017858763       -0.49997495817        0.00000013366
          6            0.50197429808       -0.16124244642       -0.00051535459        0.47170310821        0.00016564979        0.16628842255
          7            0.20164430591        0.03626697297       -0.00009425814        0.16658182358        0.00005837652       -0.47035731479
          8            0.00037645888       -0.00061933618        0.55850567097       -0.00017798963        0.50002510447        0.00000041128
          9            0.47580051857        0.14299039174       -0.17692702920       -0.23600583033        0.40815584930        0.16708658483
         10           -0.06653089742       -0.42158712194       -0.25037178283       -0.33266468146       -0.28879967263        0.23629282106
         11           -0.19738959905        0.23425408595        0.12597606171        0.28830341155        0.00011583432        0.40933305251
         12            0.47603222270        0.14259545679        0.17660774436       -0.23571450297       -0.40833028876        0.16708686476
         13           -0.06619337538       -0.42214132136        0.24952325368       -0.33287029215        0.28856503238        0.23629369780
         14            0.19755791150       -0.23453285676        0.12519167855       -0.28830326266       -0.00009062896       -0.40933374852
                          12                   13                   14
          0            0.00007554876        0.15116115423        0.26712157448
          1            0.00025227762        0.26716258890       -0.15112835540
          2           -0.30694802357        0.00025676993       -0.00005851818
          3           -0.00020381057       -0.40703205577       -0.71745691662
          4            0.00001214948        0.01272638076       -0.00717497870
          5           -0.01463937039        0.00001228602       -0.00000288351
          6            0.00019772092        0.19111368810       -0.19793342290
          7           -0.00053758796       -0.50759915745        0.58896843948
          8           -0.01465180957        0.00001228169       -0.00000304691
          9            0.22810647872       -0.15312903403       -0.00330090441
         10            0.32260980357       -0.21394411121       -0.02999429448
         11            0.54410359702       -0.39266843105       -0.03939142179
         12           -0.22836117649       -0.15274712207       -0.00338807898
         13           -0.32295520444       -0.21340401303       -0.03011786110
         14            0.54474417731        0.39175751694        0.03959935163


    $end


    """)
    testdir = 'orca_test_dir'
    file_name = 'test.hess'

    #=== Matching data values ===#
    hess = np.array([[  5.69133000e-01,  -7.00000000e-06,  -0.00000000e+00,
          -3.41891000e-01,   9.00000000e-06,   0.00000000e+00,
          -7.56840000e-02,   9.41130000e-02,   0.00000000e+00,
          -7.56850000e-02,  -4.69810000e-02,  -8.13740000e-02,
          -7.56850000e-02,  -4.69810000e-02,   8.13740000e-02],
        [ -4.00000000e-06,   5.69121000e-01,  -0.00000000e+00,
           3.00000000e-06,  -4.24950000e-02,  -0.00000000e+00,
           9.39440000e-02,  -3.08614000e-01,  -0.00000000e+00,
          -4.69710000e-02,  -1.08902000e-01,  -1.15083000e-01,
          -4.69710000e-02,  -1.08902000e-01,   1.15083000e-01],
        [ -0.00000000e+00,   0.00000000e+00,   5.69119000e-01,
           0.00000000e+00,   0.00000000e+00,  -4.24940000e-02,
          -0.00000000e+00,  -0.00000000e+00,  -4.24910000e-02,
          -8.13530000e-02,  -1.15066000e-01,  -2.41795000e-01,
           8.13530000e-02,   1.15066000e-01,  -2.41795000e-01],
        [ -3.41984000e-01,  -1.00000000e-06,   0.00000000e+00,
           3.64262000e-01,  -6.00000000e-06,   0.00000000e+00,
          -7.49300000e-03,  -1.77600000e-03,  -0.00000000e+00,
          -7.49300000e-03,   8.13000000e-04,   1.40700000e-03,
          -7.49300000e-03,   8.13000000e-04,  -1.40700000e-03],
        [  4.00000000e-06,  -4.24390000e-02,  -0.00000000e+00,
          -2.00000000e-06,   4.11320000e-02,   0.00000000e+00,
           3.07950000e-02,   2.80800000e-03,  -0.00000000e+00,
          -1.53980000e-02,  -6.49000000e-04,   2.00500000e-03,
          -1.53980000e-02,  -6.49000000e-04,  -2.00500000e-03],
        [  0.00000000e+00,  -0.00000000e+00,  -4.24330000e-02,
          -0.00000000e+00,  -0.00000000e+00,   4.11270000e-02,
           0.00000000e+00,   0.00000000e+00,  -1.84200000e-03,
          -2.66720000e-02,   2.00500000e-03,   1.66600000e-03,
           2.66720000e-02,  -2.00500000e-03,   1.66600000e-03],
        [ -7.57160000e-02,   9.41410000e-02,  -0.00000000e+00,
          -7.45800000e-03,   3.07520000e-02,   0.00000000e+00,
           7.68220000e-02,  -1.01620000e-01,  -0.00000000e+00,
           3.17800000e-03,  -1.16250000e-02,   1.42000000e-03,
           3.17800000e-03,  -1.16250000e-02,  -1.42000000e-03],
        [  9.41300000e-02,  -3.08688000e-01,   0.00000000e+00,
          -1.77800000e-03,   2.75100000e-03,  -0.00000000e+00,
          -1.01487000e-01,   3.28330000e-01,  -0.00000000e+00,
           4.57900000e-03,  -1.13170000e-02,   1.99700000e-03,
           4.57900000e-03,  -1.13180000e-02,  -1.99700000e-03],
        [ -0.00000000e+00,   0.00000000e+00,  -4.24370000e-02,
          -0.00000000e+00,  -0.00000000e+00,  -1.84100000e-03,
           0.00000000e+00,   0.00000000e+00,   4.11280000e-02,
           1.07790000e-02,  -2.44790000e-02,   1.66700000e-03,
          -1.07790000e-02,   2.44790000e-02,   1.66700000e-03],
        [ -7.57170000e-02,  -4.70670000e-02,  -8.15280000e-02,
          -7.45600000e-03,  -1.53770000e-02,  -2.66320000e-02,
           3.17700000e-03,   4.64200000e-03,   1.07500000e-02,
           7.68200000e-02,   5.07480000e-02,   8.79090000e-02,
           3.17900000e-03,   7.04600000e-03,   9.36100000e-03],
        [ -4.70650000e-02,  -1.08997000e-01,  -1.15292000e-01,
           8.89000000e-04,  -6.94000000e-04,   1.98700000e-03,
          -1.16260000e-02,  -1.12620000e-02,  -2.44470000e-02,
           5.07430000e-02,   1.12710000e-01,   1.24317000e-01,
           7.04700000e-03,   8.15900000e-03,   1.32360000e-02],
        [ -8.15190000e-02,  -1.15289000e-01,  -2.42125000e-01,
           1.53900000e-03,   1.98700000e-03,   1.60400000e-03,
           1.42200000e-03,   2.11700000e-03,   1.60300000e-03,
           8.78910000e-02,   1.24307000e-01,   2.56265000e-01,
          -9.35400000e-03,  -1.32340000e-02,  -1.78040000e-02],
        [ -7.57170000e-02,  -4.70670000e-02,   8.15290000e-02,
          -7.45600000e-03,  -1.53780000e-02,   2.66320000e-02,
           3.17800000e-03,   4.64200000e-03,  -1.07500000e-02,
           3.17900000e-03,   7.04600000e-03,  -9.36100000e-03,
           7.68200000e-02,   5.07480000e-02,  -8.79090000e-02],
        [ -4.70650000e-02,  -1.08997000e-01,   1.15292000e-01,
           8.89000000e-04,  -6.94000000e-04,  -1.98700000e-03,
          -1.16260000e-02,  -1.12620000e-02,   2.44470000e-02,
           7.04700000e-03,   8.15900000e-03,  -1.32360000e-02,
           5.07430000e-02,   1.12710000e-01,  -1.24316000e-01],
        [  8.15190000e-02,   1.15290000e-01,  -2.42125000e-01,
          -1.53900000e-03,  -1.98700000e-03,   1.60400000e-03,
          -1.42200000e-03,  -2.11700000e-03,   1.60200000e-03,
           9.35400000e-03,   1.32340000e-02,  -1.78040000e-02,
          -8.78910000e-02,  -1.24307000e-01,   2.56265000e-01]])

    atoms =['C','H','H','H','H']
    geom = np.array([-0.000000, 0.000000, 0.000000,
                            2.059801,  0.,  0.,
                            -0.6866  ,  1.942   ,  0.,
                            -0.6866  , -0.971   , -1.681821,
                            -0.6866  , -0.971   ,  1.681821
                            ])
    masses = [12.011, 1.008, 1.008, 1.008, 1.008]
    energy = -40.451555
    temp = 0.0
    freqs = np.array([    0.      ,     0.      ,     0.      ,     0.      ,
             0.      ,     0.      ,  1292.217629,  1292.583415,
          1292.939274,  1524.411513,  1524.709386,  3100.659679,
          3239.096249,  3240.029637,  3241.018763])

    modes = np.array([[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -1.23657000e-01,  -1.05760000e-02,   7.10000000e-05,
          -3.00000000e-06,   1.00000000e-06,  -1.65000000e-04,
           2.30000000e-05,   4.58120000e-02,   8.09550000e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -1.05630000e-02,   1.23653000e-01,   1.44000000e-04,
          -1.53000000e-04,  -0.00000000e+00,  -1.86000000e-04,
           7.60000000e-05,   8.09690000e-02,  -4.58020000e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -8.30000000e-05,   1.38000000e-04,  -1.24103000e-01,
          -0.00000000e+00,  -6.00000000e-06,   0.00000000e+00,
          -9.30260000e-02,   7.80000000e-05,  -1.80000000e-05],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -9.78140000e-02,  -8.36900000e-03,   5.60000000e-05,
           5.20000000e-05,   1.00000000e-06,  -4.98497000e-01,
          -2.13000000e-04,  -4.25825000e-01,  -7.50570000e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           5.13820000e-02,  -6.00708000e-01,  -7.02000000e-04,
           5.00780000e-01,   1.80000000e-04,  -1.00000000e-05,
           1.30000000e-05,   1.33140000e-02,  -7.50600000e-03],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           4.01000000e-04,  -6.72000000e-04,   6.03679000e-01,
           1.79000000e-04,  -4.99975000e-01,   0.00000000e+00,
          -1.53150000e-02,   1.30000000e-05,  -3.00000000e-06],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           5.42535000e-01,  -1.74270000e-01,  -5.57000000e-04,
           4.71703000e-01,   1.66000000e-04,   1.66288000e-01,
           2.07000000e-04,   1.99937000e-01,  -2.07069000e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           2.17937000e-01,   3.91970000e-02,  -1.02000000e-04,
           1.66582000e-01,   5.80000000e-05,  -4.70357000e-01,
          -5.62000000e-04,  -5.31035000e-01,   6.16152000e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           4.07000000e-04,  -6.69000000e-04,   6.03630000e-01,
          -1.78000000e-04,   5.00025000e-01,   0.00000000e+00,
          -1.53280000e-02,   1.30000000e-05,  -3.00000000e-06],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           5.14246000e-01,   1.54544000e-01,  -1.91222000e-01,
          -2.36006000e-01,   4.08156000e-01,   1.67087000e-01,
           2.38637000e-01,  -1.60199000e-01,  -3.45300000e-03],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -7.19070000e-02,  -4.55650000e-01,  -2.70601000e-01,
          -3.32665000e-01,  -2.88800000e-01,   2.36293000e-01,
           3.37503000e-01,  -2.23822000e-01,  -3.13790000e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -2.13339000e-01,   2.53181000e-01,   1.36154000e-01,
           2.88303000e-01,   1.16000000e-04,   4.09333000e-01,
           5.69223000e-01,  -4.10798000e-01,  -4.12090000e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           5.14496000e-01,   1.54117000e-01,   1.90877000e-01,
          -2.35715000e-01,  -4.08330000e-01,   1.67087000e-01,
          -2.38904000e-01,  -1.59800000e-01,  -3.54400000e-03],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -7.15420000e-02,  -4.56249000e-01,   2.69684000e-01,
          -3.32870000e-01,   2.88565000e-01,   2.36294000e-01,
          -3.37865000e-01,  -2.23257000e-01,  -3.15080000e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           2.13521000e-01,  -2.53482000e-01,   1.35307000e-01,
          -2.88303000e-01,  -9.10000000e-05,  -4.09334000e-01,
           5.69893000e-01,   4.09845000e-01,   4.14270000e-02]])

    dipders = np.array([[ -3.01190000e-02,  -2.50000000e-05,   0.00000000e+00],
        [  1.00000000e-05,  -3.01190000e-02,  -0.00000000e+00],
        [ -0.00000000e+00,   0.00000000e+00,  -3.01590000e-02],
        [ -1.18178000e-01,   1.90000000e-05,   0.00000000e+00],
        [ -2.90000000e-05,   7.21660000e-02,   0.00000000e+00],
        [ -0.00000000e+00,   0.00000000e+00,   7.21740000e-02],
        [  5.00740000e-02,   5.59740000e-02,  -0.00000000e+00],
        [  5.92260000e-02,  -9.72410000e-02,   0.00000000e+00],
        [  0.00000000e+00,  -0.00000000e+00,   7.21830000e-02],
        [  5.00940000e-02,  -2.79850000e-02,  -4.84610000e-02],
        [ -2.80840000e-02,   3.01860000e-02,  -6.87890000e-02],
        [ -4.87540000e-02,  -6.89430000e-02,  -4.94920000e-02],
        [  5.00940000e-02,  -2.79850000e-02,   4.84610000e-02],
        [ -2.80840000e-02,   3.01860000e-02,   6.87890000e-02],
        [  4.87530000e-02,   6.89440000e-02,  -4.94920000e-02]])

    ir_comps = np.array([[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  3.79040000e+00,   2.73700000e-01,   2.60000000e-03],
        [  3.24900000e-01,  -3.79060000e+00,  -4.30000000e-03],
        [ -2.20000000e-03,  -4.40000000e-03,   3.83080000e+00],
        [  1.29000000e-02,  -9.10000000e-03,   0.00000000e+00],
        [ -0.00000000e+00,   0.00000000e+00,   4.80000000e-03],
        [  9.15000000e-02,   1.09300000e-01,   0.00000000e+00],
        [  9.00000000e-04,   3.20000000e-03,  -3.72720000e+00],
        [  1.90370000e+00,   3.37040000e+00,   3.10000000e-03],
        [  3.50280000e+00,  -1.98250000e+00,  -7.00000000e-04]])

    ir_mags = np.array([[  0.00000000e+00],
        [  0.00000000e+00],
        [  0.00000000e+00],
        [  0.00000000e+00],
        [  0.00000000e+00],
        [  0.00000000e+00],
        [  1.44422000e+01],
        [  1.44744000e+01],
        [  1.46749000e+01],
        [  2.00000000e-04],
        [  0.00000000e+00],
        [  2.03000000e-02],
        [  1.38923000e+01],
        [  1.49837000e+01],
        [  1.61995000e+01]]).squeeze()

    polders = np.array([[ -5.94341300e+00,   2.97188100e+00,   2.97161800e+00,
          -2.47000000e-04,   0.00000000e+00,   0.00000000e+00],
        [ -1.68000000e-04,  -4.20303600e+00,   4.20264600e+00,
           2.96567900e+00,   0.00000000e+00,  -0.00000000e+00],
        [ -3.00000000e-06,  -5.00000000e-06,   2.00000000e-06,
           2.00000000e-06,   2.96532200e+00,   4.19346600e+00],
        [  9.07849800e+00,   2.17638000e-01,   2.17719000e-01,
          -5.10000000e-05,   1.00000000e-06,   0.00000000e+00],
        [  2.23000000e-04,   6.46863000e-01,  -6.46912000e-01,
           1.99675500e+00,  -0.00000000e+00,   1.00000000e-06],
        [  6.00000000e-06,  -7.00000000e-06,   1.00000000e-06,
           3.00000000e-06,   1.99668100e+00,  -6.40485000e-01],
        [ -1.04399600e+00,  -1.44161300e+00,  -6.85449000e-01,
           2.58596400e+00,   0.00000000e+00,  -1.00000000e-06],
        [  9.05935000e-01,   8.07488700e+00,  -1.13410000e-02,
          -2.03614200e+00,  -1.00000000e-06,  -2.00000000e-06],
        [ -2.30000000e-05,   1.00000000e-06,  -3.00000000e-06,
          -2.00000000e-06,  -1.26977800e+00,   1.66877300e+00],
        [ -1.04451200e+00,  -8.74492000e-01,  -1.25246500e+00,
          -1.29282500e+00,  -2.23936600e+00,  -3.30511000e-01],
        [ -4.53132000e-01,  -2.25978900e+00,  -1.77146400e+00,
          -1.46218500e+00,  -3.30621000e-01,  -2.47206900e+00],
        [ -7.84977000e-01,  -1.01697600e+00,  -5.96557900e+00,
          -3.30327000e-01,  -1.84400900e+00,  -2.60853300e+00],
        [ -1.04451700e+00,  -8.74493000e-01,  -1.25246600e+00,
          -1.29282600e+00,   2.23935900e+00,   3.30504000e-01],
        [ -4.53141000e-01,  -2.25978900e+00,  -1.77145900e+00,
          -1.46217900e+00,   3.30622000e-01,   2.47207600e+00],
        [  7.84977000e-01,   1.01696600e+00,   5.96557600e+00,
           3.30327000e-01,  -1.84401000e+00,  -2.60853100e+00]])

    raman_acts = np.array([ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
            1.5128,    1.5204,    1.5313,   26.2882,   26.2096,  140.8157,
           60.6507,   60.7749,   60.8277])

    raman_depols = np.array([0., 0.,  0.,  0.,  0.,  0.,  0.75,  0.75,  0.75,
          0.75,  0.75,  0.  ,  0.75,  0.75,  0.75])

    joblist = np.array([[ True,  True,  True],
        [ True,  True,  True],
        [ True,  True,  True],
        [ True,  True,  True],
        [ True,  True,  True]], dtype=bool)

    mwh_eigvals = np.array([ 0.     ,  0.     ,  0.     ,  0.     ,  0.     ,
          0.        ,  0.06319196,  0.06322774,  0.06326256,  0.08794171,
          0.08797609,  0.36383015,  0.39704359,  0.39727245,  0.39751505])

    mwh_eigvecs = np.array([[0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -3.94942498e-01,  -3.37783511e-02,   2.26733150e-04,
          -1.01615700e-05,   2.38076000e-06,  -5.69347300e-04,
           7.55487600e-05,   1.51161154e-01,   2.67121574e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -3.37381468e-02,   3.94929999e-01,   4.61315230e-04,
          -5.29110350e-04,  -9.81910000e-07,  -6.43013690e-04,
           2.52277620e-04,   2.67162589e-01,  -1.51128355e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -2.65229270e-04,   4.40318340e-04,  -3.96367433e-01,
          -2.16370000e-07,  -2.18289900e-05,   4.37600000e-08,
          -3.06948024e-01,   2.56769930e-04,  -5.85181800e-05],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -9.05015412e-02,  -7.74361372e-03,   5.19772600e-05,
           5.23018800e-05,   5.71510000e-07,  -4.98496537e-01,
          -2.03810570e-04,  -4.07032056e-01,  -7.17456917e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           4.75409737e-02,  -5.55800882e-01,  -6.49630810e-04,
           5.00779591e-01,   1.79653180e-04,  -9.57941000e-06,
           1.21494800e-05,   1.27263808e-02,  -7.17497870e-03],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           3.70775910e-04,  -6.21831770e-04,   5.58550833e-01,
           1.78587630e-04,  -4.99974958e-01,   1.33660000e-07,
          -1.46393704e-02,   1.22860200e-05,  -2.88351000e-06],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           5.01974298e-01,  -1.61242446e-01,  -5.15354590e-04,
           4.71703108e-01,   1.65649790e-04,   1.66288423e-01,
           1.97720920e-04,   1.91113688e-01,  -1.97933423e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           2.01644306e-01,   3.62669730e-02,  -9.42581400e-05,
           1.66581824e-01,   5.83765200e-05,  -4.70357315e-01,
          -5.37587960e-04,  -5.07599157e-01,   5.88968439e-01],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           3.76458880e-04,  -6.19336180e-04,   5.58505671e-01,
          -1.77989630e-04,   5.00025104e-01,   4.11280000e-07,
          -1.46518096e-02,   1.22816900e-05,  -3.04691000e-06],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           4.75800519e-01,   1.42990392e-01,  -1.76927029e-01,
          -2.36005830e-01,   4.08155849e-01,   1.67086585e-01,
           2.28106479e-01,  -1.53129034e-01,  -3.30090441e-03],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -6.65308974e-02,  -4.21587122e-01,  -2.50371783e-01,
          -3.32664681e-01,  -2.88799673e-01,   2.36292821e-01,
           3.22609804e-01,  -2.13944111e-01,  -2.99942945e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -1.97389599e-01,   2.34254086e-01,   1.25976062e-01,
           2.88303412e-01,   1.15834320e-04,   4.09333053e-01,
           5.44103597e-01,  -3.92668431e-01,  -3.93914218e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           4.76032223e-01,   1.42595457e-01,   1.76607744e-01,
          -2.35714503e-01,  -4.08330289e-01,   1.67086865e-01,
          -2.28361176e-01,  -1.52747122e-01,  -3.38807898e-03],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          -6.61933754e-02,  -4.22141321e-01,   2.49523254e-01,
          -3.32870292e-01,   2.88565032e-01,   2.36293698e-01,
          -3.22955204e-01,  -2.13404013e-01,  -3.01178611e-02],
        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
           1.97557911e-01,  -2.34532857e-01,   1.25191679e-01,
          -2.88303263e-01,  -9.06289600e-05,  -4.09333749e-01,
           5.44744177e-01,   3.91757517e-01,   3.95993516e-02]])



    #=== Defining ways to break the .hess ===#
    class names(object):
        hess = 'hess'
        geom = 'geom'
        atsym = 'atsym'
        energy = 'energy'
        temp = 'temp'
        freqs = 'freq'
        modes = 'modes'
        dipders = 'dipders'
        ir = 'ir'
        polders = 'polders'
        raman = 'raman'
        joblist = 'joblist'
        mwh_eigvals = 'mwh_eigvals'
        mwh_eigvecs = 'mwh_eigvecs'

        suffix_dim2 = '_dim2'
        suffix_badfreq = '_badfreq'
        suffix_badval = '_badval'

        E = frozenset([hess, geom, atsym, energy, temp, freqs, modes,
                    dipders, ir, polders, raman, mwh_eigvals, mwh_eigvecs])
        suffixes = frozenset([suffix_dim2, suffix_badfreq, suffix_badval])


    bad_block_substs = {
            names.hess  : ('$hessian', '$harshman'),
            names.geom  : ('$atoms', '$martians'),
            names.energy : ('$act_energy', '$fact_harbaly'),
            names.temp  : ('$actual_temp', '$schmactuish'),
            names.freqs : ('$vibrational_freq', '$varbifishing_sweg'),
            names.modes : ('$normal_modes', '$formal_toads'),
            names.dipders : ('$dipole_deriv', '$tadpole_gurriv'),
            names.ir    : ('$ir_spectrum', '$lidar_fulcrum'),
            names.polders : ('$polarizability_derivatives',
                                            '$polarbear_skivs'),
            names.raman : ('$raman_spectrum', '$jammin_sparkrum'),
            names.joblist : ('$job_list', '$mob_gist'),
            names.mwh_eigvals : ('$eigenvalues_mass_weighted_hessian',
                                        '$goy_that_is_long'),
            names.mwh_eigvecs : ('$eigenvectors_mass_weighted_hessian',
                                        '$also_very_long')
                        }

    trunc_block_substs = {
            names.hess  : ('7       0.094130  -0.308688', 'gabrab'),
            names.hess + names.suffix_dim2 :
                            ('10       0.007047   0.008159', 'farfrif'),
            names.geom  : ('H      1.0080      2.059801', 'gaffraf'),
            names.freqs : ('10     1524.709386', 'fobbardgik'),
            names.modes : ('1      -0.010563   0.123653', 'gommerbik'),
            names.dipders: ('0.050074     0.055974', 'cheezenugget'),
            names.ir    : ('1292.94      14.6749', 'spuffle'),
            names.polders: ('9.078498     0.217638', 'frabbitz'),
            names.raman : ('1524.71      26.2096', 'zorbotz'),
            names.joblist : ('3    1    1    1', 'gorb'),
            names.mwh_eigvals: ('8          0.06326255701', 'zabbert'),
            names.mwh_eigvecs: ('9            0.47580051857', 'fobborj')
                            }

    bad_data_substs = {
            names.hess  : ('$hessian\n15', '$hessian\n30'),
            names.freqs : ('al_frequencies\n15', 'al_frequencies\n30'),
            names.atsym : ('C     12.0110     -0.000000',
                                    'Cx    12.0110     -0.000000'),
            names.modes : ('l_modes\n15', 'l_modes\n19'),
            names.modes + names.suffix_dim2 :
                            ('l_modes\n15 15', 'l_modes\n15 19'),
            names.dipders : ('ole_derivatives\n15', 'ole_derivatives\n25'),
            names.dipders + names.suffix_badval :
                            ('-0.118178     0.000019', '-11817.8     0.000019'),
            names.ir    : ('ir_spectrum\n15', 'ir_spectrum\n38'),
            names.ir + names.suffix_badfreq :
                            ('1292.94      14.6749', '3232.28      14.6749'),
            names.polders : ('ity_derivatives\n15', 'ity_derivatives\n38'),
            names.raman : ('man_spectrum\n15', 'man_spectrum\n27'),
            names.raman + names.suffix_badfreq :
                            ('1524.71      26.2096', '3184.71      26.2096'),
            names.joblist : ('job_list\n15', 'job_list\n38'),
            names.mwh_eigvals: ('values_mass_weighted_hessian\n15',
                                        'values_mass_weighted_hessian\n7'),
            names.mwh_eigvecs: ('ectors_mass_weighted_hessian\n15',
                                        'ectors_mass_weighted_hessian\n39'),
            names.mwh_eigvecs + names.suffix_dim2:
                                ('ectors_mass_weighted_hessian\n15 15',
                                        'ectors_mass_weighted_hessian\n15 19')
                        }

    alt_data_substs = {
            names.joblist : ('2    1    1    1', '2    0    1    0')
                        }

    alt_data_values = {
            names.joblist : np.matrix([[ True,  True,  True],
                                        [ True,  True,  True],
                                        [ False,  True,  False],
                                        [ True,  True,  True],
                                        [ True,  True,  True]], dtype=bool)
                        }

## end class SuperORCAHess


class TestORCAHessKnownGood(SuperORCAHess):
    # Testing errors related to file parsing of a known-good HESS file

    @classmethod
    def setUpClass(self):

        setUpTestDir(self.testdir)

        # Write the file
        with open(self.file_name, 'w') as f:
            f.write(self.file_text_good)


    @classmethod
    def tearDownClass(self):
        import os

        # Delete the hess file
        os.remove(self.file_name)

        # Remove the test directory
        tearDownTestDir(self.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.hess import ORCA_HESS

        # Create the object
        self.oh = ORCA_HESS(self.file_name)

        # Enable long messages
        self.longMessage = True

    def test_HESS_KnownGoodAtomVec(self):
        self.assertEqual(len(self.oh.atom_syms), len(self.atoms))
        for i in range(len(self.oh.atom_syms)):
            self.assertEqual(self.oh.atom_syms[i], self.atoms[i])

    def test_HESS_KnownGoodAtomMasses(self):
        self.assertEqual(len(self.oh.atom_masses), len(self.masses))
        for i in range(len(self.oh.atom_masses)):
            self.assertAlmostEqual(self.oh.atom_masses[i],
                        self.masses[i], delta=1e-4,
                        msg="Atom index " + str(i))

    def test_HESS_KnownGoodHess(self):
        self.assertEqual(self.oh.hess.shape, self.hess.shape)
        for i in range(self.oh.hess.shape[0]):
            for j in range(self.oh.hess.shape[1]):
                self.assertAlmostEqual(self.oh.hess[i,j],
                            self.hess[i,j], delta=1e-6,
                            msg="Hessian element (" + str(i) + ',' +
                                                        str(j) + ')')

    def test_HESS_KnownGoodGeom(self):
        self.assertEqual(self.oh.geom.shape, self.geom.shape)
        for i in range(self.oh.geom.shape[0]):
            self.assertAlmostEqual(self.oh.geom[i],
                        self.geom[i], delta=1e-6,
                        msg="Coordinate index " + str(i))

    def test_HESS_KnownGoodCheckGeomWorks(self):
        self.assertTrue(self.oh.check_geom(self.geom, self.atoms))

    def test_HESS_KnownGoodEnergy(self):
        self.assertAlmostEqual(self.oh.energy, self.energy, delta=1e-6)

    def test_HESS_KnownGoodTemp(self):
        self.assertAlmostEqual(self.oh.temp, self.temp, delta=1e-6)

    def test_HESS_KnownGoodFreqs(self):
        self.assertEqual(self.oh.freqs.shape, self.freqs.shape)
        for i in range(self.oh.freqs.shape[0]):
            self.assertAlmostEqual(self.oh.freqs[i],
                        self.freqs[i], delta=1e-6,
                        msg="Frequency index " + str(i))

    def test_HESS_KnownGoodInitFlagDefined(self):
        self.assertTrue('initialized' in self.oh.__dict__)
        if 'initialized' in self.oh.__dict__:
            self.assertTrue(self.oh.initialized)

    def test_HESS_KnownGoodModes(self):
        self.assertEqual(self.oh.modes.shape, self.modes.shape)
        for i in range(self.oh.modes.shape[0]):
            for j in range(self.oh.modes.shape[1]):
                self.assertAlmostEqual(self.oh.modes[i,j],
                            self.modes[i,j], delta=1e-6,
                            msg="Mode " + str(j) + ", element " + str(i))

    def test_HESS_KnownGoodDipDers(self):
        self.assertEqual(self.oh.dipders.shape, self.dipders.shape)
        for i in range(self.oh.dipders.shape[0]):
            for j in range(self.oh.dipders.shape[1]):
                self.assertAlmostEqual(self.oh.dipders[i,j],
                            self.dipders[i,j], delta=1e-6,
                            msg="Dipole derivative element (" + str(i) + ',' +
                                                        str(j) + ')')

    def test_HESS_KnownGoodIRSpectrum(self):
        self.assertEqual(self.oh.ir_comps.shape, self.ir_comps.shape)
        self.assertEqual(self.oh.ir_mags.shape, self.ir_mags.shape)
        for i in range(self.oh.ir_comps.shape[0]):
            self.assertAlmostEqual(self.oh.ir_mags[i], self.ir_mags[i],
                    delta=1e-4, msg="IR T**2 element (" + str(i) + ')')
            for j in range(self.oh.ir_comps.shape[1]):
                self.assertAlmostEqual(self.oh.ir_comps[i,j],
                            self.ir_comps[i,j], delta=1e-4,
                            msg="IR T_i derivative element (" + str(i) + ',' +
                                                        str(j) + ')')

    def test_HESS_KnownGoodPolDers(self):
        self.assertEqual(self.oh.polders.shape, self.polders.shape)
        for i in range(self.oh.polders.shape[0]):
            for j in range(self.oh.polders.shape[1]):
                self.assertAlmostEqual(self.oh.polders[i,j],
                            self.polders[i,j], delta=1e-4,
                            msg="Polarizability derivative element (" +
                                        str(i) + ',' + str(j) + ')')

    def test_HESS_KnownGoodRamanSpectrum(self):
        self.assertEqual(self.oh.raman_acts.shape, self.raman_acts.shape)
        self.assertEqual(self.oh.raman_depols.shape, self.raman_depols.shape)
        for i in range(self.oh.raman_acts.shape[0]):
            self.assertAlmostEqual(self.oh.raman_acts[i],
                    self.raman_acts[i],
                    delta=1e-4, msg="Raman activity element (" + str(i) + ')')
            self.assertAlmostEqual(self.oh.raman_depols[i],
                        self.raman_depols[i], delta=1e-4,
                        msg="Raman depolarization element (" + str(i) + ')')

    def test_HESS_KnownGoodJobList(self):
        self.assertEqual(self.oh.joblist.shape, self.joblist.shape)
        for i in range(self.oh.joblist.shape[0]):
            for j in range(3):
                self.assertEqual(self.oh.joblist[i,j], self.joblist[i,j],
                        msg="Job list element (" + str(i) + ',' + str(j) + ')')

    def test_HESS_KnownGoodMWHEigvals(self):
        self.assertEqual(self.oh.mwh_eigvals.shape, self.mwh_eigvals.shape)
        for i in range(self.oh.mwh_eigvals.shape[0]):
            self.assertAlmostEqual(self.oh.mwh_eigvals[i],
                    self.mwh_eigvals[i],
                    delta=1e-8, msg="MWH eigenvector (" + str(i) + ')')

    def test_HESS_KnownGoodMWHEigvecs(self):
        self.assertEqual(self.oh.mwh_eigvecs.shape, self.mwh_eigvecs.shape)
        for i in range(self.oh.mwh_eigvecs.shape[0]):
            for j in range(self.oh.mwh_eigvecs.shape[1]):
                self.assertAlmostEqual(self.oh.mwh_eigvecs[i,j],
                        self.oh.mwh_eigvecs[i,j],
                        delta=1e-8, msg="MWH eigenvectors element (" +
                                                str(i) + ',' + str(j) + ')')


## end class TestORCAEngradKnownGood


class TestORCAHessMissingBlocks(SuperORCAHess):
    # Ensuring importing a non-HESS file throws the right errors


    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.bad_block_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.bad_block_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any .hess files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in
                                            self.bad_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_HESS_MissingBlockHess(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.hess_block, self.file_name +
                    self.names.hess)

    def test_HESS_MissingBlockGeom(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.at_block, self.file_name +
                    self.names.geom)

    def test_HESS_MissingBlockEnergy(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.energy, self.file_name +
                    self.names.energy)

    def test_HESS_MissingBlockTemp(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.temp, self.file_name +
                    self.names.temp)

    def test_HESS_MissingBlockFreqs(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.freq_block, self.file_name +
                    self.names.freqs)

    def test_HESS_MissingBlockModes(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.modes_block, self.file_name +
                    self.names.modes)

    def test_HESS_MissingBlockDipders(self):

        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name +
                                            self.names.dipders).dipders)

    def test_HESS_MissingBlockIRSpectrum(self):

        from opan.hess import ORCA_HESS

        h = ORCA_HESS(self.file_name + self.names.ir)

        self.assertIsNone(h.ir_comps)
        self.assertIsNone(h.ir_mags)

    def test_HESS_MissingBlockPolders(self):

        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name +
                                            self.names.polders).polders)

    def test_HESS_MissingBlockRamanSpectrum(self):

        from opan.hess import ORCA_HESS

        h = ORCA_HESS(self.file_name + self.names.raman)

        self.assertIsNone(h.raman_acts)
        self.assertIsNone(h.raman_depols)

    def test_HESS_MissingBlockJobList(self):

        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name +
                                            self.names.joblist).joblist)

    def test_HESS_MissingBlockMWHEigvals(self):

        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name +
                                        self.names.mwh_eigvals).mwh_eigvals)

    def test_HESS_MissingBlockMWHEigvecs(self):

        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name +
                                        self.names.mwh_eigvecs).mwh_eigvecs)

## end class TestORCAHessMissingBlocks


class TestORCAHessTruncatedBlocks(SuperORCAHess):
    # Ensuring importing a HESS file with an incomplete block throws the
    #       right errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.trunc_block_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.trunc_block_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in
                                            self.trunc_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_HESS_TruncatedBlocksHess(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # 'Early' (non-final section) truncation
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.hess_block, self.file_name + self.names.hess)
        # 'Late' (final section) truncation
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.hess_block,
                    self.file_name + self.names.hess + self.names.suffix_dim2)

    def test_HESS_TruncatedBlocksGeom(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.at_block, self.file_name + self.names.geom)

    def test_HESS_TruncatedBlocksFreqs(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.freq_block, self.file_name + self.names.freqs)

    def test_HESS_TruncatedBlocksModes(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.modes_block, self.file_name + self.names.modes)

    def test_HESS_TruncatedBlocksDipders(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.dipder_block, self.file_name + self.names.dipders)

    def test_HESS_TruncatedBlocksIRSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.ir_block, self.file_name + self.names.ir)

    def test_HESS_TruncatedBlocksPolders(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.polder_block, self.file_name + self.names.polders)

    def test_HESS_TruncatedBlocksRamanSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.raman_block, self.file_name + self.names.raman)

    def test_HESS_TruncatedBlocksJobList(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.job_block, self.file_name + self.names.joblist)

    def test_HESS_TruncatedBlocksMWHEigvals(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.eigval_block,
                    self.file_name + self.names.mwh_eigvals)

    def test_HESS_TruncatedBlocksMWHEigvecs(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.eigvec_block,
                    self.file_name + self.names.mwh_eigvecs)

## end class TestORCAHessTruncatedBlocks


class TestORCAHessBadData(SuperORCAHess):
    # Ensuring importing a HESS file with data of valid formatting but
    #  invalid content raises the appropriate errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.bad_data_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.bad_data_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + dname) for dname in
                                            self.bad_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_HESS_BadDataHessDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.hess_block, self.file_name + self.names.hess)

    def test_HESS_BadDataFreqDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.freq_block, self.file_name + self.names.freqs)

    def test_HESS_BadDataAtomSym(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        self.assertRaises(KeyError, ORCA_HESS,
                                    self.file_name + self.names.atsym)

    def test_HESS_BadDataNormalModes(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # First dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.modes_block, self.file_name + self.names.modes)
        # Second dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.modes_block,
                    self.file_name + self.names.modes + self.names.suffix_dim2)

    def test_HESS_BadDataDipdersDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.dipder_block, self.file_name + self.names.dipders)

    def test_HESS_BadDataDipdersTooBig(self):

        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name + self.names.dipders +
                        self.names.suffix_badval).dipders)

    def test_HESS_BadDataIRSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # Bad dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.ir_block, self.file_name + self.names.ir)
        # Mismatched frequency
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.ir_block,
                    self.file_name + self.names.ir + self.names.suffix_badfreq)

    def test_HESS_BadDataPoldersDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.polder_block, self.file_name + self.names.polders)

    def test_HESS_BadDataRamanSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # Bad dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.raman_block, self.file_name + self.names.raman)
        # Mismatched frequency
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.raman_block,
                    self.file_name + self.names.raman +
                                                self.names.suffix_badfreq)

    def test_HESS_BadDataJobListDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.job_block, self.file_name + self.names.joblist)

    def test_HESS_BadDataMWHEigvalsDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.eigval_block,
                    self.file_name + self.names.mwh_eigvals)

    def test_HESS_BadDataMWHEigvecs(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # First dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.eigvec_block,
                    self.file_name + self.names.mwh_eigvecs)
        # Second dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS,
                    HESSError.eigvec_block,
                    self.file_name + self.names.mwh_eigvecs +
                            self.names.suffix_dim2)

## end class TestORCAHessBadData


class TestORCAHessAltData(SuperORCAHess):
    # Ensuring importing a HESS file with data of valid formatting but
    #  invalid content raises the appropriate errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.alt_data_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.alt_data_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + dname) for dname in
                                            self.alt_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def setUp(self):
        # Long messages
        self.longMessage = True

    def test_HESS_AltDataJobList(self):

        from opan.hess import ORCA_HESS

        h = ORCA_HESS(self.file_name + self.names.joblist)

        self.assertEqual(h.joblist.shape,
                        self.alt_data_values[self.names.joblist].shape)
        for i in range(h.joblist.shape[0]):
            for j in range(3):
                self.assertEqual(h.joblist[i,j],
                        self.alt_data_values[self.names.joblist][i,j],
                        msg="Job list element (" + str(i) + ',' + str(j) + ')')

## end class TestORCAHessAltData


class TestORCAHessLiveData(SuperORCAHess):

    def setUp(self):
        # Enable long messages
        self.longMessage = True

    def test_HESS_LiveData(self):

        import os
        from opan.hess import ORCA_HESS
        from opan.error import HESSError

        for fname in os.listdir(self.resourcedir):
            if fname[:9] == "test_orca" and fname[-4:] == "hess":
                print("\nTesting file '" + fname + "' ....")
                try:
                    ORCA_HESS(os.path.join(self.resourcedir, fname))
                except (IOError, HESSError) as e: # pragma: no cover
                    self.fail("Load of test file '" + str(fname) +
                            "' failed:\n" + str(e))

## end class TestORCAHessLiveData


# =============================  OPAN_XYZ  ================================== #

class SuperOPANXYZ(unittest.TestCase):
    # Superclass for all OPAN .xyz test cases

    # Imports
    from textwrap import dedent
    import numpy as np

    # Superclass constants

    file_text_good = dedent("""\
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.753188968196
      Cu     -0.698514     -0.206760     -0.206944
      O       1.004955      0.146285      0.146415
      H       1.346594      1.037721     -0.077335
      H       1.346963     -0.077245      1.037865
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.755551440250
      Cu     -0.718110     -0.216851     -0.216979
      O       0.992982      0.155412      0.155540
      H       1.362425      1.052057     -0.090702
      H       1.362702     -0.090619      1.052142
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.757560388801
      Cu     -0.749991     -0.233409     -0.233458
      O       0.976820      0.174023      0.174112
      H       1.386513      1.058149     -0.098827
      H       1.386658     -0.098763      1.058173
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.758833365916
      Cu     -0.782499     -0.250038     -0.250011
      O       0.963852      0.194315      0.194346
      H       1.409311      1.054489     -0.098804
      H       1.409336     -0.098766      1.054470
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759333332418
      Cu     -0.803394     -0.259878     -0.259811
      O       0.960032      0.206083      0.206067
      H       1.421697      1.047777     -0.094003
      H       1.421665     -0.093982      1.047748
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759472543226
      Cu     -0.812328     -0.262882     -0.262805
      O       0.962913      0.207978      0.207941
      H       1.424728      1.045012     -0.090123
      H       1.424687     -0.090109      1.044988
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759565737953
      Cu     -0.817928     -0.263432     -0.263355
      O       0.969218      0.205351      0.205307
      H       1.424373      1.045141     -0.087073
      H       1.424337     -0.087060      1.045122
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759656625423
      Cu     -0.822757     -0.262377     -0.262309
      O       0.979838      0.198680      0.198643
      H       1.421470      1.048067     -0.084385
      H       1.421448     -0.084370      1.048052
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759687014064
      Cu     -0.823630     -0.261011     -0.260952
      O       0.985847      0.194139      0.194114
      H       1.418895      1.051207     -0.084355
      H       1.418887     -0.084335      1.051194
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691012126
      Cu     -0.823025     -0.260523     -0.260468
      O       0.986742      0.193258      0.193240
      H       1.418142      1.052463     -0.085220
      H       1.418140     -0.085197      1.052449
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691138943
      Cu     -0.822700     -0.260487     -0.260435
      O       0.986486      0.193462      0.193448
      H       1.418106      1.052644     -0.085643
      H       1.418108     -0.085620      1.052631
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691121397
      Cu     -0.822625     -0.260504     -0.260454
      O       0.986384      0.193616      0.193603
      H       1.418117      1.052652     -0.085788
      H       1.418122     -0.085764      1.052640
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691135231
      Cu     -0.822616     -0.260519     -0.260472
      O       0.986370      0.193707      0.193696
      H       1.418117      1.052648     -0.085861
      H       1.418127     -0.085836      1.052638
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691151054
      Cu     -0.822627     -0.260520     -0.260475
      O       0.986393      0.193712      0.193703
      H       1.418111      1.052648     -0.085866
      H       1.418123     -0.085840      1.052639
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691154491
      Cu     -0.822630     -0.260516     -0.260472
      O       0.986410      0.193700      0.193692
      H       1.418103      1.052651     -0.085861
      H       1.418117     -0.085835      1.052642
    4
    Coordinates from ORCA-job CuII_1aq_128 E -1715.759691151236
      Cu     -0.822629     -0.260515     -0.260471
      O       0.986409      0.193698      0.193691
      H       1.418102      1.052651     -0.085860
      H       1.418117     -0.085834      1.052642
    """)

    testdir = 'opan_test_dir'
    file_name = 'test.xyz'

    class names(object):
        num_atoms = 'num_atoms'
        num_geoms = 'num_geoms'
        file_start = 'file_start'
        const_numats = 'const_numats'
        atom_counts = 'atom_counts'

        bad_lategeom_atomnum = 'bad_lategeom_atomnum'
        diff_lategeom_atomnum = 'diff_lategeom_atomnum'
        bad_1stgeom_atomnum = 'bad_1stgeom_atomnum'

        bad_lategeom_atomsym = 'bad_lategeom_atomsym'
        diff_lategeom_atomsym = 'diff_lategeom_atomsym'
        bad_1stgeom_atomsym = 'bad_1stgeom_atomsym'

        Cu_1st_as_atomnum = 'Cu_1st_as_atomnum'
        Cu_late_as_atomnum = 'Cu_late_as_atomnum'
    ## end class names

    # Known-good data for this particular trajectory
    num_atoms = 4
    num_geoms = 16

    atom_syms = ['CU', 'O', 'H', 'H']
    geoms = [np.array([-1.32000015, -0.39071977, -0.39106748,  1.89908972,  0.27643859,
        0.27668425,  2.54469386,  1.96100848, -0.14614197,  2.54539117,
       -0.14597189,  1.9612806 ]),
         np.array([-1.35703123, -0.409789  , -0.41003088,  1.87646403,  0.29368612,
                0.293928  ,  2.57461012,  1.9880996 , -0.17140194,  2.57513357,
               -0.17124509,  1.98826022]),
         np.array([-1.41727759, -0.44107909, -0.44117168,  1.84592227,  0.32885581,
                0.329024  ,  2.62012984,  1.99961181, -0.18675596,  2.62040385,
               -0.18663502,  1.99965716]),
         np.array([-1.4787088 , -0.47250334, -0.47245232,  1.8214163 ,  0.36720213,
                0.36726071,  2.66321181,  1.99269541, -0.1867125 ,  2.66325906,
               -0.18664069,  1.99265951]),
         np.array([-1.51819463, -0.49109825, -0.49097163,  1.81419755,  0.38944043,
                0.38941019,  2.68661796,  1.98001157, -0.17763992,  2.68655749,
               -0.17760024,  1.97995677]),
         np.array([-1.53507744, -0.49677498, -0.49662947,  1.81964185,  0.39302146,
                0.39295154,  2.69234572,  1.97478648, -0.17030779,  2.69226824,
               -0.17028133,  1.97474112]),
         np.array([-1.54565991, -0.49781433, -0.49766882,  1.83155657,  0.38805715,
                0.387974  ,  2.69167487,  1.97503025, -0.16454412,  2.69160684,
               -0.16451956,  1.97499435]),
         np.array([-1.5547854 , -0.49582067, -0.49569217,  1.85162547,  0.37545079,
                0.37538087,  2.68618899,  1.98055959, -0.15946454,  2.68614742,
               -0.15943619,  1.98053124]),
         np.array([-1.55643513, -0.49323931, -0.49312781,  1.86298083,  0.36686954,
                0.3668223 ,  2.68132295,  1.98649333, -0.15940785,  2.68130783,
               -0.15937005,  1.98646876]),
         np.array([-1.55529184, -0.49231712, -0.49221318,  1.86467214,  0.36520469,
                0.36517068,  2.67989999,  1.98886683, -0.16104246,  2.67989621,
               -0.160999  ,  1.98884037]),
         np.array([-1.55467768, -0.49224909, -0.49215082,  1.86418837,  0.3655902 ,
                0.36556374,  2.67983196,  1.98920887, -0.16184181,  2.67983574,
               -0.16179835,  1.9891843 ]),
         np.array([-1.55453595, -0.49228121, -0.49218673,  1.86399561,  0.36588121,
                0.36585665,  2.67985274,  1.98922398, -0.16211582,  2.67986219,
               -0.16207047,  1.98920131]),
         np.array([-1.55451895, -0.49230956, -0.49222074,  1.86396916,  0.36605318,
                0.36603239,  2.67985274,  1.98921643, -0.16225377,  2.67987164,
               -0.16220653,  1.98919753]),
         np.array([-1.55453973, -0.49231145, -0.49222641,  1.86401262,  0.36606263,
                0.36604562,  2.6798414 ,  1.98921643, -0.16226322,  2.67986408,
               -0.16221409,  1.98919942]),
         np.array([-1.5545454 , -0.49230389, -0.49222074,  1.86404475,  0.36603995,
                0.36602483,  2.67982629,  1.98922209, -0.16225377,  2.67985274,
               -0.16220464,  1.98920509]),
         np.array([-1.55454351, -0.492302  , -0.49221885,  1.86404286,  0.36603617,
                0.36602294,  2.6798244 ,  1.98922209, -0.16225189,  2.67985274,
               -0.16220275,  1.98920509])]

    dist_Cu_O = np.array([3.354629, 3.383183, 3.440127, 3.507285, 3.557423,
                3.582909, 3.602044, 3.622363, 3.629285, 3.628582, 3.627700,
                3.627540, 3.627595, 3.627663, 3.627685, 3.627680])
    dist_O_H1 = np.array([1.852933, 1.890762, 1.912290, 1.912522, 1.900678,
                1.892315, 1.887740, 1.886512, 1.889388, 1.891501, 1.891975,
                1.891988, 1.891935, 1.891909, 1.891903, 1.891904])

    angle_Cu_O_H1 = np.array([118.017, 119.211, 121.278, 123.250, 124.283,
                124.437, 124.193, 123.537, 123.041, 122.911, 122.909, 122.911,
                122.910, 122.909, 122.907, 122.907])

    dihed_H2_O_Cu_H1 = np.array([131.287, 135.462, 142.180, 148.963, 152.551,
                152.735, 151.323, 148.425, 146.671, 146.400, 146.499, 146.558,
                146.589, 146.589, 146.584, 146.584])


    good_direct_geom = [np.array([-1.32000015, -0.39071977, -0.39106748,
         1.89908972,  0.27643859, 0.27668425,  2.54469386,  1.96100848,
         -0.14614197,  2.54539117, -0.14597189,  1.9612806 ])]
    good_direct_atoms = ['CU', 'O', 'H', 'H']
    good_direct_O1Dist = np.array([3.3546284826543507, 0.0,
                                1.8529334749323307, 1.8531057866658869])
    good_direct_CuO1Angle = np.array([0.0, None, 118.0167027446415,
                                            118.05344488961961], dtype=object)
    good_direct_Dihed = np.array([None, None, None, 131.28736349459098],
                                                                dtype=object)

    bad_file_data_substs = {
                names.file_start : ('4\nCoordin', '\n4\nCoordin'),
                names.atom_counts :
                            ('H       1.421470      1.048067', 'shmarb'),
                names.const_numats :
                            ('1.052449\n4', '1.052449\n6'),
                names.bad_lategeom_atomnum :
                            ('Cu     -0.822629', '300    -0.822629'),
                names.diff_lategeom_atomnum :
                            ('Cu     -0.822629', '45     -0.822629'),
                names.bad_1stgeom_atomnum :
                            ('Cu     -0.698514', '300     -0.698514'),
                names.bad_lategeom_atomsym :
                            ('Cu     -0.822629', 'Xd     -0.822629'),
                names.diff_lategeom_atomsym :
                            ('Cu     -0.822629', 'Se     -0.822629'),
                names.bad_1stgeom_atomsym :
                            ('Cu     -0.698514', 'Xd      -0.698514')
                            }

    alt_file_data_substs = {
                names.Cu_1st_as_atomnum :
                            ('Cu     -0.698514', '29     -0.822629'),
                names.Cu_late_as_atomnum :
                            ('Cu     -0.822629', '29     -0.822629')
                            }

## end class SuperOPANXYZ


class TestOPANXYZGoodFileData(SuperOPANXYZ):
    # Ensuring importing a known OpenBabel xyz file with good data reports
    #  the correct geometric parameters, etc.

    @classmethod
    def setUpClass(self):
        # Set up the directory and add the good file
        setUpTestDir(self.testdir)

        # Write the file
        with open(self.file_name, 'w') as f:
            f.write(self.file_text_good)

    @classmethod
    def tearDownClass(self):
        import os

        # Delete the xyz file
        os.remove(self.file_name)

        # Remove the test directory
        tearDownTestDir(self.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.xyz import OPAN_XYZ

        # Create the object
        self.xyz = OPAN_XYZ(path=self.file_name)

        # Long messages
        self.longMessage = True

    def test_XYZ_GoodFileDataNumAtoms(self):
        self.assertEqual(self.xyz.num_atoms, self.num_atoms)

    def test_XYZ_GoodFileDataNumGeoms(self):
        self.assertEqual(self.xyz.num_geoms, self.num_geoms)

    def test_XYZ_GoodFileDataAtomSyms(self):
        for i in range(len(self.xyz.atom_syms)):
            self.assertEqual(self.xyz.atom_syms[i],
                    self.atom_syms[i],
                    msg="Coordinates element (" + str(i) + ')')

    def test_XYZ_GoodFileDataCoords(self):
        for g in range(len(self.geoms)):
            for i in range(self.geoms[g].shape[0]):
                self.assertAlmostEqual(self.xyz.geoms[g][i],
                        self.geoms[g][i],
                        delta=1e-6,
                        msg="Geometry #" + str(g) +
                                ", coordinate element #" + str(i))

    def test_XYZ_GoodFileDataDistances(self):
        for t in zip(self.dist_Cu_O, self.xyz.Dist_iter(None, 0, 1),
                                        range(self.dist_Cu_O.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-5,
                        msg="Cu-O distance mismatch at geom #" + str(t[2]) +
                        ": " + str(t[0:2]))
        for t in zip(self.dist_O_H1, self.xyz.Dist_iter(None, 1, 2),
                                        range(self.dist_O_H1.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-5,
                        msg="O-H1 distance mismatch at geom #" + str(t[2]) +
                        ": " + str(t[0:2]))

    def test_XYZ_GoodFileDataAngles(self):
        for t in zip(self.angle_Cu_O_H1, self.xyz.Angle_iter(None, 0, 1, 2),
                                        range(self.angle_Cu_O_H1.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-2,
                        msg="Cu-O-H1 angle mismatch at geom #" + str(t[2]) +
                        ": " + str(t[0:2]))

    def test_XYZ_GoodFileDataDihedrals(self):
        for t in zip(self.dihed_H2_O_Cu_H1,
                                    self.xyz.Dihed_iter(None, 3, 1, 0, 2),
                                    range(self.dihed_H2_O_Cu_H1.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-2,
                        msg="H2-Cu-O-H1 dihedral mismatch at geom #" +
                        str(t[2]) + ": " + str(t[0:2]))

    def test_XYZ_GoodFileDataIterGeom(self):

        import numpy as np

        idxs = [1,4,8]
        for t in zip(np.array(self.geoms)[idxs],
                        self.xyz.Geom_iter(idxs),
                        range(len(idxs))
                    ):
            for i in range(t[0].shape[0]):
                self.assertAlmostEqual(t[0][i], t[1][i], delta=1e-6,
                        msg="Geometry #" + str(t[2]) + \
                                ", coordinate element #" + str(i))

    #TEST: Displ_iter call (Dist_iter only calls Displ_single)

## end class TestOPANXYZGoodData


class TestOPANXYZAltFileData(SuperOPANXYZ):
    # Ensuring successful import of an XYZ file with valid data of alternative
    #  formatting to that contained in SuperOPANXYZ.file_text_good
    #  (primarily if elements are specified by atomic number).

    @classmethod
    def setUpClass(self):
        # Set up the directory and add alternate-g files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.alt_file_data_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.alt_file_data_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any created files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in
                                            self.alt_file_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_XYZ_AltFileDataFirstCuAtomnum(self):

        from opan.xyz import OPAN_XYZ
        from opan.error import XYZError

        try:
            xyz = OPAN_XYZ(path=(self.file_name + self.names.Cu_1st_as_atomnum))
        except XYZError:  # pragma: no cover
            self.fail("XYZ import failed when success was expected.")

        self.assertEqual(xyz.atom_syms[0].upper(),
                                        self.good_direct_atoms[0])

    def test_XYZ_AltFileDataLateCuAtomnum(self):

        from opan.xyz import OPAN_XYZ
        from opan.error import XYZError

        try:
            xyz = OPAN_XYZ(path=(self.file_name +
                                            self.names.Cu_late_as_atomnum))
        except XYZError:  # pragma: no cover
            self.fail("XYZ import failed when success was expected.")

        self.assertEqual(xyz.atom_syms[0].upper(),
                                        self.good_direct_atoms[0])

## end class TestOPANXYZAltFileData


class TestOPANXYZGoodDirectData(SuperOPANXYZ):
    # Confirming sanity of an OPAN_XYZ generated directly from data.

    def setUp(self):
        # Load the object

        # Imports
        from opan.xyz import OPAN_XYZ

        # Create the object
        self.xyz = OPAN_XYZ(atom_syms=self.good_direct_atoms,
                                            coords=self.good_direct_geom)

        # Long messages
        self.longMessage = True

    def test_XYZ_GoodDirectDataNumAtoms(self):
        self.assertEqual(self.xyz.num_atoms, len(self.good_direct_atoms))

    def test_XYZ_GoodDirectDataNumGeoms(self):
        self.assertEqual(self.xyz.num_geoms, 1)

    def test_XYZ_GoodDirectDataNumCoords(self):
        self.assertEqual(self.xyz.geoms[0].shape[0], self.xyz.num_atoms * 3)

    def test_XYZ_GoodDirectDataCoords(self):
        for i in range(self.good_direct_geom[0].shape[0]):
            self.assertAlmostEqual(self.xyz.geoms[0][i],
                    self.good_direct_geom[0][i],
                    delta=1e-8,
                    msg="Coordinates element (" + str(i) + ')')

    def test_XYZ_GoodDirectDataIterO1Dist(self):
        for tup in zip(
                    self.good_direct_O1Dist,
                    self.xyz.Dist_iter(0,1,None),
                    range(self.good_direct_O1Dist.shape[0])
                        ):
            self.assertAlmostEqual(tup[0], tup[1], delta=1e-5,
                    msg="Distance between O1 and atom #" + str(tup[2]) +
                            " (" + self.xyz.atom_syms[tup[2]].capitalize() +
                            ")")

    def test_XYZ_GoodDirectDataIterCuO1Angle(self):
        for tup in zip(
                    self.good_direct_CuO1Angle,
                    self.xyz.Angle_iter(0,0,1,None),
                    range(self.good_direct_CuO1Angle.shape[0])
                        ):
            self.assertAlmostEqual(tup[0], tup[1], delta=5e-3,
                    msg="Angle Cu-O-X with atom #" + str(tup[2]) +
                            " (" + self.xyz.atom_syms[tup[2]].capitalize() +
                            ")")

    def test_XYZ_GoodDirectDataIterDihed(self):
        for tup in zip(
                    self.good_direct_Dihed,
                    self.xyz.Dihed_iter(0,None, 1,0,2),
                    range(self.good_direct_Dihed.shape[0])
                        ):
            self.assertAlmostEqual(tup[0], tup[1], delta=5e-3,
                    msg="Dihedral with atom #" + str(tup[2]) +
                            " (" + self.xyz.atom_syms[tup[2]].capitalize() +
                            ")")

## end class TestOPANXYZGoodDirectData


class TestOPANXYZBadFileData(SuperOPANXYZ):
    # Ensuring importing an XYZ file with data of generally valid formatting but
    #  invalid content raises the appropriate errors

    @classmethod
    def setUpClass(self):
        # Set up the directory and add munged files

        setUpTestDir(self.testdir)

        # Write the files
        for bname in self.bad_file_data_substs.keys():
            with open(self.file_name + bname, 'w') as f:
                f.write(self.file_text_good
                                    .replace(*self.bad_file_data_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any created files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in
                                            self.bad_file_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_XYZ_BadFileDataFileStart(self):
        # File without a number-of-atoms spec at the very start should throw
        #  an XYZError

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.file_start))

    def test_XYZ_BadFileDataAtomCounts(self):
        # Multiple-geometry file without a consistent number of atoms
        #  throughout should throw an XYZError, regardless of the 'numats'
        #  spec at the top of each geometry block

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.atom_counts))

    def test_XYZ_BadFileDataConstNumats(self):
        # Multiple-geometry file should have identical num_ats specs at the
        #  head of each geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.const_numats))

    def test_XYZ_BadFileDataBadLateGeomAtomNum(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.bad_lategeom_atomnum))

    def test_XYZ_BadFileDataDifferentLateGeomAtomNum(self):
        # Discrepant atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                    XYZError.xyzfile,
                    path=(self.file_name + self.names.diff_lategeom_atomnum))

    def test_XYZ_BadFileDataBadFirstGeomAtomNum(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.bad_1stgeom_atomnum))

    def test_XYZ_BadFileDataBadLateGeomAtomSym(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.bad_lategeom_atomsym))

    def test_XYZ_BadFileDataDifferentLateGeomAtomSym(self):
        # Discrepant atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                    XYZError.xyzfile,
                    path=(self.file_name + self.names.diff_lategeom_atomsym))

    def test_XYZ_BadFileDataBadFirstGeomAtomSym(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        assertErrorAndTypecode(self, XYZError, OPAN_XYZ,
                        XYZError.xyzfile,
                        path=(self.file_name + self.names.bad_1stgeom_atomsym))

## end class TestOPANXYZBadFileData


class TestOPANXYZBadDirectData(SuperOPANXYZ):
    # Ensuring invalid data passed directly to construct an OPAN_XYZ instance
    #  throws the appropriate errors or otherwise exhibits the correct problems.

    def test_XYZ_BadDirectDataTruncCoords(self):
        from opan.xyz import OPAN_XYZ
        self.assertRaises(ValueError, OPAN_XYZ, \
                    coords=self.good_direct_geom[0][:-2], \
                    atom_syms=self.good_direct_atoms)

    def test_XYZ_BadDirectDataBadElement(self):
        from opan.xyz import OPAN_XYZ

        # Copy symbols, munge, and pass into assert
        munge_atoms = self.good_direct_atoms[:]
        munge_atoms[0] = 'CX'
        self.assertRaises(ValueError, OPAN_XYZ, \
                    coords=self.good_direct_geom, \
                    atom_syms=munge_atoms)

    def test_XYZ_BadDirectDataNonRealCoords(self):
        from opan.xyz import OPAN_XYZ
        import numpy as np

        # Copy coords, munge, and pass to assert
        munge_coords = self.good_direct_geom[0].copy() * 1.j
        self.assertRaises(ValueError, OPAN_XYZ, \
                    coords=munge_coords, \
                    atom_syms=self.good_direct_atoms)

    # TEST: ValueError for Angle at_1 == at_2 and at_3 == at_2 cases

    # TEST: ValueError for Dihed any at_x equal

    # TEST: IndexError for invalid (out-of-range) at_x for Dist, Angle, Dihed

    # TEST: XYZError.dihed for too-nearly-linear atom trio(s)

    # TEST: ValueError when multiple 'None' values passed to an X_iter method

    # TEST: ValueError when 'None' passed to an X_iter with another iterable

## end class TestOPANXYZBadDirectData


class TestOPANXYZBadUsage(SuperOPANXYZ):
    # Ensuring bad use cases of OPAN_XYZ throw the appropriate errors or
    #  otherwise exhibit the correct problems.

    @classmethod
    def setUpClass(self):
        # Set up the directory and add the good file
        setUpTestDir(self.testdir)

        # Write the file
        with open(self.file_name, 'w') as f:
            f.write(self.file_text_good)

    @classmethod
    def tearDownClass(self):
        import os

        # Delete the xyz file
        os.remove(self.file_name)

        # Remove the test directory
        tearDownTestDir(self.testdir)

    def test_XYZ_BadUsageBadInitParams(self):
        # Calling without 'path' or both 'atom_syms' and 'coords'

        from opan.xyz import OPAN_XYZ

        self.assertRaises(NameError, OPAN_XYZ)
        self.assertRaises(NameError, OPAN_XYZ, atom_syms=['CU', 'O'])
        self.assertRaises(NameError, OPAN_XYZ, coords=[1,2,3])

    def test_XYZ_BadUsageReInitDirectData(self):
        # Should pitch a fit if someone uses the private data method to try to
        #  reinit an OPAN_XYZ instance

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        # Should work fine
        xyz = OPAN_XYZ(coords=self.good_direct_geom, \
                                atom_syms=self.good_direct_atoms)

        # Should raise error to re-init with same contents.
        assertErrorAndTypecode(self, XYZError, xyz._load_data, \
                        XYZError.overwrite, \
                        coords=self.good_direct_geom, \
                        atom_syms=self.good_direct_atoms)

    def test_XYZ_BadUsageReInitFileData(self):
        # Should pitch a fit if someone uses the private file method to try to
        #  reinit an OPAN_XYZ instance

        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        # Should work fine
        xyz = OPAN_XYZ(path=self.file_name)

        # Should raise error to re-init with same contents.
        assertErrorAndTypecode(self, XYZError, xyz._load_file, \
                        XYZError.overwrite, \
                        XYZ_path=self.file_name)

    def test_XYZ_BadUsageNotXYZFile(self):

        import os
        from opan.error import XYZError
        from opan.xyz import OPAN_XYZ

        badfilename = 'bad.file'

        try:
            with open(badfilename, 'w') as f:
                f.write("This is not an OpenBabel file.\n\n")
                f.write("(In case you were wondering...)\n\n")

            assertErrorAndTypecode(self, XYZError, OPAN_XYZ, \
                            XYZError.xyzfile, \
                            path=badfilename)

        finally:
            os.remove(badfilename)


## end class TestOPANXYZBadUsage


# ===========================  utils.inertia  =============================== #

class SuperOPANUtilsInertia(object):
    # Superclass for all test cases for the various top types. Implemented
    #  as a top-level class rather than a subclass of unittest.TestCase to
    #  prevent the tests being run spuriously on the superclass. Also avoids
    #  run-time errors since self.fname does not exist in the superclass
    #  namespace.
    #
    # Subclasses must also subclass unittest.TestCase, and must declare
    #  SuperOPANUtilsInertia as a higher priority superclass. Otherwise, the
    #  default form of the `setUpClass` method from unittest.TestCase will be
    #  called.

    # Imports
    import os

    # Common constants
    filedir = os.path.join('resource','inertia')

    # Common setUpClass method
    @classmethod
    def setUpClass(self):

        # Imports
        import os
        from opan.xyz import OPAN_XYZ as XYZ
        from opan.hess import ORCA_HESS as HESS

        # Attach class objects
        basename = os.path.join(self.filedir, self.fname)
        self.xyz = XYZ(path=basename + '.xyz')
        self.hess = HESS(basename + '.hess')

        # Always long messages
        self.longMessage = True

    def test_ctr_mass(self):
        import opan.utils.inertia as oui
        ctr_mass = oui.ctr_mass(self.xyz.geoms[0], self.hess.atom_masses)
        for i in range(3):
            self.assertAlmostEqual(self.ctr_mass[i], ctr_mass[i],
                        delta=1e-7,
                        msg="Center-of-mass index '" + str(i) + "'")

    def test_ctr_geom(self):
        import opan.utils.inertia as oui
        ctr_geom = oui.ctr_geom(self.xyz.geoms[0], self.hess.atom_masses)
        for i in range(ctr_geom.shape[0]):
            self.assertAlmostEqual(self.ctr_geom[i], ctr_geom[i],
                        delta=1e-7,
                        msg="Centered geometry index '" + str(i) + "'")

    def test_i_tensor(self):
        import opan.utils.inertia as oui
        i_tensor = oui.inertia_tensor(self.xyz.geoms[0],
                                                    self.hess.atom_masses)
        for i in range(i_tensor.shape[0]):
            for j in range(i_tensor.shape[1]):
                self.assertAlmostEqual(self.i_tensor[i,j],
                            i_tensor[i,j],
                            delta=1e-7,
                            msg="Inertia tensor element (" + str(i) + "," +
                                                             str(j) + ")")

    def test_moments(self):
        import opan.utils.inertia as oui
        moments = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[0]
        for i in range(moments.shape[0]):
            if moments[i] < 1e-10:
                self.assertAlmostEqual(self.moments[i], moments[i],
                        delta=1e-7,
                        msg="Principal moment index '" + str(i) + "'")
            else:
                self.assertAlmostEqual(self.moments[i] / moments[i], 1.0,
                        delta=1e-7,
                        msg="Principal moment index '" + str(i) + "'")

    def test_axes(self):
        import opan.utils.inertia as oui
        axes = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[1]
        for i in range(axes.shape[0]):
            for j in range(axes.shape[1]):
                if abs(axes[i,j]) < 1e-10:
                    self.assertAlmostEqual(self.axes[i,j], axes[i,j],
                            delta=1e-7,
                            msg="Principal axis #" + str(j) + ", element " +
                                                                        str(i))
                else:
                    self.assertAlmostEqual(self.axes[i,j] / axes[i,j], 1.0,
                            delta=1e-7,
                            msg="Principal axis #" + str(j) + ", element " +
                                                                        str(i))

    def test_axes_orthonorm(self):
        import opan.utils.inertia as oui
        from opan.utils.vector import orthonorm_check as onchk
        axes = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[1]
        on, nfail, ofail = onchk(axes, report=True)
        self.assertTrue(on, msg="Norm failures: " + str(nfail) +
                                "; ortho failures: " + str(ofail))

    def test_rot_consts(self):
        import opan.utils.inertia as oui
        from opan.const import PRM, EU_RotConst as EURC
        rc = oui.rot_consts(self.xyz.geoms[0], self.hess.atom_masses,
                                                            EURC.InvInertia)
        for i in range(rc.shape[0]):
            if rc[i] >= 1/(2.0*PRM.Zero_Moment_Tol):
                self.assertAlmostEqual(self.rc[i], rc[i],
                        delta=1.0,
                        msg="Rotational constant index '" + str(i) + "'")
            else:
                self.assertAlmostEqual(self.rc[i] / rc[i], 1.0,
                        delta=PRM.Equal_Moment_Tol,
                        msg="Rotational constant index '" + str(i) + "'")

    def test_toptype(self):
        import opan.utils.inertia as oui
        top = oui.principals(self.xyz.geoms[0], self.hess.atom_masses)[2]
        self.assertEqual(top, self.top)

## end class SuperOPANUtilsInertia


class TestOPANUtilsInertiaAsymm(SuperOPANUtilsInertia, unittest.TestCase):
    # Asymmetric molecule test-case (H2O). Includes proofreading tests for,
    #  e.g., bad object shapes

    # Imports
    import numpy as np
    from opan.const import E_TopType, EU_RotConst as EURC

    # Constants for superclass method use
    fname = "H2O_Asymm"
    ctr_mass = np.array([-11.42671795,   0.2152039 ,  -0.02145012])
    ctr_geom = np.array([-0.07503316, -0.09586698,  0.04549049,
                1.7685783 ,  0.01137086, -0.00539566,
                -0.57765012,  1.51023212, -0.71663048])
    i_tensor = np.array([[ 2.99702139,  0.74400967, -0.35304507],
                [ 0.74400967,  4.13012039,  1.16077065],
                [-0.35304507,  1.16077065,  6.02553148]])
    moments = np.array([ 2.41469379,  4.16164284,  6.57633663])
    axes = np.array([[  8.16492638e-01,  -5.77355846e-01,   6.07765301e-15],
                [ -5.21610191e-01,  -7.37657520e-01,  -4.28700585e-01],
                [  2.47512789e-01,   3.50030872e-01,  -9.03446627e-01]])
    top = E_TopType.Asymmetrical
    rc = np.array([ 0.207065592447,  0.120144860861,  0.076030171231])

    rc_units = {EURC.AngFreqAtomic: np.array([  1.135920230159e-04,   6.590905634729e-05,   4.170862410476e-05]),
         EURC.AngFreqSeconds: np.array([  4.696050231566e+12,   2.724770904719e+12,   1.724291800473e+12]),
         EURC.CyclicFreqAtomic: np.array([  1.807873195879e-05,   1.048975211219e-05,   6.638133695834e-06]),
         EURC.CyclicFreqHz: np.array([  7.473996073616e+11,   4.336607582790e+11,   2.744295633781e+11]),
         EURC.CyclicFreqMHz: np.array([ 747399.60736156872 ,  433660.758279046626,  274429.563378054823]),
         EURC.InvInertia: np.array([ 0.207065592447,  0.120144860861,  0.076030171231]),
         EURC.WaveNumAtomic: np.array([  1.318826101078e-07,   7.652173233680e-08,   4.842454659134e-08]),
         EURC.WaveNumCM: np.array([ 24.922201369646,  14.460511669382,   9.150913076387])}


    def test_UtilsInertiaBadGeomShape(self):
        from opan.utils.inertia import ctr_mass as cm, _fadnOv as fO, \
                                                        _fadnPv as fP
        import numpy as np
        self.assertRaises(ValueError, cm, np.array([[1,1],[1,1]]),
                                                    self.hess.atom_masses)
        self.assertRaises(ValueError, fO, [1,2,3], self.xyz.geoms[0][:-2])
        self.assertRaises(ValueError, fP, [1,2,3], self.xyz.geoms[0][:-2])

    def test_UtilsInertiaBadMassesShape(self):
        from opan.utils.inertia import ctr_mass as cm
        import numpy as np
        self.assertRaises(ValueError, cm, self.xyz.geoms[0],
                                                     np.array([[1,1],[1,1]]))

    def test_UtilsInertiaBadGeomLength(self):
        from opan.utils.inertia import ctr_mass as cm
        import numpy as np
        self.assertRaises(ValueError, cm, self.xyz.geoms[0][:-2],
                                                    self.hess.atom_masses)

    def test_UtilsInertiaBadGeomMassesLengths(self):
        from opan.utils.inertia import ctr_mass as cm
        import numpy as np
        self.assertRaises(ValueError, cm, self.xyz.geoms[0],
                                                self.hess.atom_masses[:-2])

    def test_UtilsInertiaBadRefVecShape(self):
        from opan.utils.inertia import _fadnOv as fO, _fadnPv as fP
        import numpy as np
        self.assertRaises(ValueError, fO, [1,2,3,4], self.xyz.geoms[0])
        self.assertRaises(ValueError, fP, [1,2,3,4], self.xyz.geoms[0])

    def test_UtilsInertiaBadRefVecNorm(self):
        from opan.utils.inertia import _fadnOv as fO, _fadnPv as fP
        import numpy as np
        self.assertRaises(ValueError, fO, [1e-7, 0, 1e-9], self.xyz.geoms[0])
        self.assertRaises(ValueError, fP, [1e-7, 0, 1e-9], self.xyz.geoms[0])

    def test_UtilsInertiaCheckAllRotConstUnits(self):
        from opan.utils.inertia import rot_consts
        from opan.const import EU_RotConst as EURC
        for (u, a) in [(u, rot_consts(self.xyz.geoms[0],
                    self.hess.atom_masses, units=u)) for u in EURC]:
            for i in range(a.shape[0]):
                self.assertAlmostEqual(self.rc_units[u][i] / a[i], 1.0,
                        delta=1e-8,
                        msg="Rotational constant units '" + str(u) +
                                                ",' index '" + str(i) + "'")

    def test_UtilsInertiaBadRotConstUnits(self):
        from opan.utils.inertia import rot_consts
        self.assertRaises(ValueError, rot_consts, self.xyz.geoms[0],
                                self.hess.atom_masses, units="ThisIsInvalid")

## end class TestOPANUtilsInertiaAsymm


class TestOPANUtilsInertiaAtom(SuperOPANUtilsInertia, unittest.TestCase):
    # Lone atom test-case (Cu).

    # Imports
    import numpy as np
    from opan.const import E_TopType, PRM

    # Constants for superclass method use
    fname = "Cu_Atom"
    ctr_mass = np.array([-1.88972612,  3.77945225, -7.5589045 ])
    ctr_geom = np.array([0.0, 0.0, 0.0])
    i_tensor = np.zeros((3,3))
    moments = np.zeros((3,))
    axes = np.eye(3)
    top = E_TopType.Atom
    rc = np.repeat(1.0/(2.0*PRM.Zero_Moment_Tol), 3)

## end class TestOPANUtilsInertiaAtom


class TestOPANUtilsInertiaLinear(SuperOPANUtilsInertia, unittest.TestCase):
    # Linear molecule test-case (chloroethyne).

    # Imports
    import numpy as np
    from opan.const import E_TopType, PRM

    # Constants for superclass method use
    fname = "HC2Cl_Linear"
    ctr_mass = np.array([ 3.78157026, -0.13914054,  0.11918503])
    ctr_geom = np.array([ -1.81226617e+00,   4.52442528e-10,   7.10376657e-09,
                         1.29696963e+00,   4.73253989e-09,  -1.02091483e-08,
                         3.58244725e+00,  -6.42275907e-09,  -1.13833863e-08,
                         5.59880599e+00,   4.22696256e-09,   7.43858825e-09])
    i_tensor = np.array([[  5.44288079e-15,   2.07854690e-07,   1.06328882e-06],
                       [  2.07854690e-07,   3.22388376e+02,  -4.43485080e-16],
                       [  1.06328882e-06,  -4.43485080e-16,   3.22388376e+02]])
    moments = np.array([  5.68434189e-14,   3.22388376e+02,   3.22388376e+02])
    axes = np.array([[ -1.00000000e+00,   2.49655671e-10,  -3.91982518e-09],
                   [  2.49655671e-10,   1.00000000e+00,   0.00000000e+00],
                   [  3.91982518e-09,  -9.78606588e-19,  -1.00000000e+00]])
    top = E_TopType.Linear
    rc = np.array([1.0/(2.0*PRM.Equal_Moment_Tol),
                                    1.550924404326e-03, 1.550924404326e-03])

    def test_UtilsInertiaLinearNoNonParallelVec(self):
        from opan.utils.inertia import _fadnPv as fP, ctr_geom as cg
        import numpy as np
        from opan.error import INERTIAError as INErr
        assertErrorAndTypecode(self, INErr, fP , INErr.bad_geom,
                            self.xyz.Displ_single(0,0,1),
                            cg(self.xyz.geoms[0], self.hess.atom_masses))

## end class TestOPANUtilsInertiaLinear


class TestOPANUtilsInertiaSymmProl(SuperOPANUtilsInertia, unittest.TestCase):
    # Prolate symmetric test case (chloromethane)

    # Imports
    import numpy as np
    from opan.const import E_TopType

    # Constants for superclass method use
    fname = "CH3Cl_SymmProl"
    ctr_mass = np.array([ 1.19526288,  1.19526288,  1.19526288])
    ctr_geom = np.array([-1.36880665, -1.36880665, -1.36880665,
                            0.61213832,  0.61213832, 0.61213832,
                            -2.54039904, -2.54039904, -0.13884799,
                            -0.13884799, -2.54039904, -2.54039904,
                            -2.54039904, -0.13884799, -2.54039905])
    i_tensor = np.array([[ 97.63769388, -43.0052599 , -43.00525991],
                       [-43.0052599 ,  97.63769391, -43.00525989],
                       [-43.00525991, -43.00525989,  97.63769389]])
    moments = np.array([  11.6271741 ,  140.64295379,  140.6429538 ])
    axes = np.array([[ -5.77350269e-01,  -4.08248290e-01,  -7.07106781e-01],
                       [ -5.77350269e-01,  -4.08248291e-01,   7.07106781e-01],
                       [ -5.77350269e-01,   8.16496581e-01,   2.78286699e-10]])
    top = E_TopType.SymmProlate
    rc = np.array([ 0.043002710352,  0.003555101671,  0.003555101671])

## end class TestOPANUtilsInertiaSymmProl


class TestOPANUtilsInertiaSymmObl(SuperOPANUtilsInertia, unittest.TestCase):
    # Oblate symmetric test case (ammonia)

    # Imports
    import numpy as np
    from opan.const import E_TopType

    # Constants for superclass method use
    fname = "NH3_SymmObl"
    ctr_mass = np.array([ -1.36098776e+01,   1.70206822e+00,  -7.62922335e-04])
    ctr_geom = np.array([ -4.79117194e-02,  -1.35356483e-01,   5.88277405e-03,
                         1.89623378e+00,   3.55505121e-02,  -1.54206207e-03,
                        -6.15221552e-01,   8.55891028e-01,  -1.57668545e+00,
                        -6.15238957e-01,   9.89449583e-01,   1.49648146e+00])
    i_tensor = np.array([[ 6.74683794,  0.98560398, -0.04281628],
                           [ 0.98560398,  9.1833784 , -0.12106317],
                           [-0.04281628, -0.12106317,  6.40284824]])
    moments = np.array([ 6.39758672,  6.39807449,  9.53740337])
    axes = np.array([[-0.941357532609, -0.052449650272, -0.333309210773],
                   [ 0.334894464181, -0.024822665809, -0.94192862422 ],
                   [ 0.041130203771, -0.998315015137,  0.040932100963]])
    top = E_TopType.SymmOblate
    rc = np.array([ 0.078154470132,  0.078148511783,  0.05242517071 ])

## end class TestOPANUtilsInertiaSymmObl


class TestOPANUtilsInertiaPlanar(SuperOPANUtilsInertia, unittest.TestCase):
    # Planar (oblate symmetric special case) test case (benzene)

    # Imports
    import numpy as np
    from opan.const import E_TopType

    # Constants for superclass method use
    fname = "C6H6_Planar"
    ctr_mass = np.array(
            [0.0, 0.0, 0.0])
    ctr_geom = np.array([ 1.889726124565,  0.            ,  0.,
                0.944863062283, 1.636550830098,  0.,
                -0.944863062283,  1.636550830098, 0.,
                -1.889726124565,  0.            ,  0.,
                -0.944863062283, -1.636550830098,  0.,
                0.944863062283, -1.636550830098,  0.,
                3.77945224913 ,  0., 0.,
                1.889726124565,  3.273101660137,  0.,
                -1.889726124565,  3.273101660137,  0.,
                -3.77945224913 , 0.,  0.,
                -1.889726124565, -3.273101660137, 0.,
                 1.889726124565, -3.273101660137,  0.])
    i_tensor = np.array([[ 171.871779008624,   -0.,   -0.],
                            [  -0.,  171.871779003993,   -0.],
                            [  -0., -0.,  343.743558012617]])
    moments = np.array([171.871779003993, 171.871779008624, 343.743558012617])
    axes = np.array([[ 1.,  0.,  0.],
                       [ 0.,  1.,  0.],
                       [ 0.,  0.,  1.]])
    top = E_TopType.SymmOblate
    rc = np.array([ 0.002909145427,  0.002909145427,  0.001454572714])

## end class TestOPANUtilsInertiaPlanar


class TestOPANUtilsInertiaSpher(SuperOPANUtilsInertia, unittest.TestCase):
    # Spherical test case (methane)

    # Imports
    import numpy as np
    from opan.const import E_TopType

    # Constants for superclass method use
    fname = "CH4_Spher"
    ctr_mass = np.array(
            [  8.734135537996e-11,   4.387156003093e-11,  -1.188484647673e-11])
    ctr_geom = np.array(
            [  4.000728815448e-11,   2.009566928560e-11,  -5.443942085537e-12,
             1.180503717314e+00,   1.180503717408e+00,   1.180503717381e+00,
            -1.180503717657e+00,  -1.180503717520e+00,   1.180503717479e+00,
             1.180503717435e+00,  -1.180503717544e+00,  -1.180503717454e+00,
            -1.180503717568e+00,   1.180503717417e+00,  -1.180503717341e+00])
    i_tensor = np.array(
            [[  1.123790191290e+01,   7.775846633251e-11,   5.015445836420e-10],
           [  7.775846633251e-11,   1.123790191310e+01,  -3.777778090353e-11],
           [  5.015445836420e-10,  -3.777778090353e-11,   1.123790191366e+01]])
    moments = np.array([ 11.237901912633,  11.237901913117,  11.237901913909])
    axes = np.array(
            [[  5.773502691632e-01,  -4.082482905123e-01,   7.071067811802e-01],
           [  5.773502692095e-01,  -4.082482904247e-01,  -7.071067811929e-01],
           [  5.773502691962e-01,   8.164965809231e-01,   6.946399011554e-11]])
    top = E_TopType.Spherical
    rc = np.array([ 0.044492290811,  0.044492290809,  0.044492290806])

## end class TestOPANUtilsInertiaSpher


# ==========================  Helper Functions  ============================= #

def assertErrorAndTypecode(testclass, errtype, cobj, tc, *args, **kwargs):
    """ Wrapper for asserting correct OPANErrors and proper typecodes.

    Function tests (using testclass.assertX methods) whether 'cobj' raises
    'errtype' with typecode 'tc' when instantiated/called with *args and
    **kwargs.

    Parameters
    ----------
    testclass   : object reference
        Subclass of unittest.TestCase (or related), from which the .assertX
        methods should be called
    errtype     : object reference
        Subclass of OPANError expected to be raised
    cobj        : object reference
        Callable object to be instantiated or called
    tc          : str / typecode "enum"
        Typecode to check for
    *args and **kwargs are passed to the instantiation of 'cobj'

    Returns
    -------
    (none)

    """

    # Assert the proper error
    testclass.assertRaises(errtype, cobj, *args, **kwargs)

    # Ensure correct typecode; suppress repeated error, and ignore any
    #  other error raised, as it will have been reported by the above
    #  .assertRaises call
    try:
        out = cobj(*args, **kwargs)
    except errtype as err:
        testclass.assertEqual(err.tc, tc)
    except:  # pragma: no cover
        pass

## end def assertErrorAndTypecode


def setUpTestDir(dirname):
    """ Create and change working directory to test directory.

    Parameters
    ----------
    dirname : str
        Name of desired working directory
    """

    import os

    # Check if test directory already exists (or file of same name);
    #  error if so
    if os.path.isdir(dirname) or os.path.isfile(dirname): # pragma: no cover
        raise(IOError("Cannot create new test directory!"))

    # Create and change to test directory
    os.mkdir(dirname)
    os.chdir(dirname)


def tearDownTestDir(dirname):
    """ Exit and attempt removal of test directory

    Parameters
    ----------
    dirname: str
        Name of working directory
    """

    import os

    # Switch to parent directory
    os.chdir(os.path.pardir)

    # Try to remove the temp directory
    os.rmdir(dirname)


if __name__ == '__main__':

    # Run test package, default verbose. Can override with '-q' at
    #  commandline
    unittest.main(verbosity=2)

## end main block
