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
    # Superclass for all ORCA engrad test cases

    # Imports
    from textwrap import dedent
    import numpy as np

    # Superclass constants
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

    atoms = np.array([['CU'],['O'],['H'],['H']])
    gradient = np.matrix([0.000004637000, 0.000001922807, 0.000002366827, \
                0.000010704036, -0.000004735732, -0.000005580943, \
                0.000011997130, -0.000008391237, -0.000007912155, \
                0.000011493781, -0.000008177479, -0.000008233658]).transpose()
    energy = -1715.759691151236
    geom = np.matrix([-1.5545432, -0.4923021, -0.4922193, \
                            1.8640436, 0.3660366, 0.3660223, \
                            2.6798241, 1.9892213, -0.1622520, \
                            2.6798521, -0.1622023, 1.9892043]).transpose()

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

    def test_ENGRAD_KnownGoodAtomVec(self):
        # Confirm the values coming out of the ENGRAD match the known-good
        #  example file.

        # Confirm the atom list is good
        for i in range(self.oe.atom_syms.shape[0]):
            self.assertEqual(self.oe.atom_syms[i,0], self.atoms[i,0])

    def test_ENGRAD_KnownGoodGradient(self):
        # Confirm the known-good gradient matches what's expected.

        # Confirm the gradient vector is good
        self.longMessage = True
        for i in range(self.oe.gradient.shape[0]):
            self.assertAlmostEqual(self.oe.gradient[i,0], \
                        self.gradient[i,0], delta=1e-10, \
                        msg="Gradient index " + str(i))

    def test_ENGRAD_KnownGoodEnergy(self):
        # Confirm the known-good energy matches what's expected
        self.assertAlmostEqual(self.energy, self.oe.energy, delta=1e-12)

    def test_ENGRAD_KnownGoodGeom(self):
        self.longMessage = True
        for i in range(self.oe.geom_vec.shape[0]):
            self.assertAlmostEqual(self.oe.geom_vec[i,0], \
                        self.geom[i,0], delta=1e-7, \
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
                f.write(self.file_text_good \
                                    .replace(*self.bad_block_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in \
                                            self.bad_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_ENGRAD_MissingBlockNumAts(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.numats, self.file_name + self.names.numats)

    def test_ENGRAD_MissingBlockEnergy(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.en, self.file_name + self.names.energy)

    def test_ENGRAD_MissingBlockGrad(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.gradblock, self.file_name + self.names.grad)

    def test_ENGRAD_MissingBlockGeom(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
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
                f.write(self.file_text_good \
                                    .replace(*self.trunc_block_substs[bname]))


    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in \
                                        self.trunc_block_substs.keys()]

        # Try to remove the directory
        tearDownTestDir(self.testdir)

    def test_ENGRAD_TruncatedBlockGrad(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.gradblock, self.file_name + self.names.grad)

    def test_ENGRAD_TruncatedBlockGeom(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
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
                f.write(self.file_text_good \
                                    .replace(*self.bad_data_substs[dname]))

    @classmethod
    def tearDownClass(self):
        # Tear down test directory and remove files

        import os

        # Try to remove the files
        [os.remove(self.file_name + dname) for dname in \
                                            self.bad_data_substs.keys()]

        # Remove the test directory
        tearDownTestDir(self.testdir)


    def test_ENGRAD_BadDataAtomicNum(self):
        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.geomblock, self.file_name \
                                        + self.names.atomicnum)

    def test_ENGRAD_BadDataNumAtoms(self):
        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.gradblock, self.file_name \
                                        + self.names.numats)

## end class TestORCAEngradBadData



# ============================  ORCA_ERROR ================================== #

class TestOPANErrorInitErrors(unittest.TestCase):
    # Testing errors that should be thrown on initialization

    def test_OPANError_init_NotImplemented(self):
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

## end class TestXYZErrorInitConfig


# =============================  ORCA_HESS  ================================= #

class SuperORCAHess(unittest.TestCase):
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
    hess = np.matrix([[  5.69133000e-01,  -7.00000000e-06,  -0.00000000e+00,
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

    atoms = np.array([['C'],['H'],['H'],['H'],['H']])
    geom = np.matrix([-0.000000, 0.000000, 0.000000, \
                            2.059801,  0.,  0., \
                            -0.6866  ,  1.942   ,  0., \
                            -0.6866  , -0.971   , -1.681821, \
                            -0.6866  , -0.971   ,  1.681821
                            ]).transpose()
    masses = np.matrix([12.011, 1.008, 1.008, 1.008, 1.008]).transpose()
    energy = -40.451555
    temp = 0.0
    freqs = np.matrix([[    0.      ,     0.      ,     0.      ,     0.      ,
             0.      ,     0.      ,  1292.217629,  1292.583415,
          1292.939274,  1524.411513,  1524.709386,  3100.659679,
          3239.096249,  3240.029637,  3241.018763]]).transpose()

    modes = np.matrix([[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
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

    dipders = np.matrix([[ -3.01190000e-02,  -2.50000000e-05,   0.00000000e+00],
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

    ir_comps = np.matrix([[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
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

    ir_mags = np.matrix([[  0.00000000e+00],
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
        [  1.61995000e+01]])

    polders = np.matrix([[ -5.94341300e+00,   2.97188100e+00,   2.97161800e+00,
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

        suffix_dim2 = '_dim2'
        suffix_badfreq = '_badfreq'

        E = frozenset([hess, geom, atsym, energy, temp, freqs, modes, \
                    dipders, ir, polders])
        suffixes = frozenset([suffix_dim2, suffix_badfreq])


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
                                            '$polarbear_skivs')
                        }

    trunc_block_substs = {
            names.hess  : ('7       0.094130  -0.308688', 'gabrab'),
            names.geom  : ('H      1.0080      2.059801', 'gaffraf'),
            names.freqs : ('10     1524.709386', 'fobbardgik'),
            names.modes : ('1      -0.010563   0.123653', 'gommerbik'),
            names.dipders: ('0.050074     0.055974', 'cheezenugget'),
            names.ir    : ('1292.94      14.6749', 'spuffle'),
            names.polders: ('9.078498     0.217638', 'frabbitz')
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
            names.ir    : ('ir_spectrum\n15', 'ir_spectrum\n38'),
            names.ir + names.suffix_badfreq :
                            ('1292.94      14.6749', '3232.28      14.6749'),
            names.polders : ('ity_derivatives\n15', 'ity_derivatives\n38')
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

    def test_HESS_KnownGoodAtomVec(self):
        self.assertEqual(self.oh.atom_syms.shape, self.atoms.shape)
        for i in range(self.oh.atom_syms.shape[0]):
            self.assertEqual(self.oh.atom_syms[i,0], self.atoms[i,0])

    def test_HESS_KnownGoodAtomMasses(self):
        self.longMessage = True
        self.assertEqual(self.oh.atom_masses.shape, self.masses.shape)
        for i in range(self.oh.atom_masses.shape[0]):
            self.assertAlmostEqual(self.oh.atom_masses[i,0], \
                        self.masses[i,0], delta=1e-4, \
                        msg="Atom index " + str(i))

    def test_HESS_KnownGoodHess(self):
        self.longMessage = True
        self.assertEqual(self.oh.hess.shape, self.hess.shape)
        for i in range(self.oh.hess.shape[0]):
            for j in range(self.oh.hess.shape[1]):
                self.assertAlmostEqual(self.oh.hess[i,j], \
                            self.hess[i,j], delta=1e-6, \
                            msg="Hessian element (" + str(i) + ',' + \
                                                        str(j) + ')')

    def test_HESS_KnownGoodGeom(self):
        self.longMessage = True
        self.assertEqual(self.oh.geom.shape, self.geom.shape)
        for i in range(self.oh.geom.shape[0]):
            self.assertAlmostEqual(self.oh.geom[i,0], \
                        self.geom[i,0], delta=1e-6, \
                        msg="Coordinate index " + str(i))

    def test_HESS_KnownGoodEnergy(self):
        self.assertAlmostEqual(self.oh.energy, self.energy, delta=1e-6)

    def test_HESS_KnownGoodTemp(self):
        self.assertAlmostEqual(self.oh.temp, self.temp, delta=1e-6)

    def test_HESS_KnownGoodFreqs(self):
        self.longMessage = True
        self.assertEqual(self.oh.freqs.shape, self.freqs.shape)
        for i in range(self.oh.freqs.shape[0]):
            self.assertAlmostEqual(self.oh.freqs[i,0], \
                        self.freqs[i,0], delta=1e-6, \
                        msg="Frequency index " + str(i))

    def test_HESS_KnownGoodInitFlagDefined(self):
        self.assertTrue('initialized' in self.oh.__dict__)
        if 'initialized' in self.oh.__dict__:
            self.assertTrue(self.oh.initialized)

    def test_HESS_KnownGoodModes(self):
        self.longMessage = True
        self.assertEqual(self.oh.modes.shape, self.modes.shape)
        for i in range(self.oh.modes.shape[0]):
            for j in range(self.oh.modes.shape[1]):
                self.assertAlmostEqual(self.oh.modes[i,j], \
                            self.modes[i,j], delta=1e-6, \
                            msg="Mode " + str(j) + ", element " + str(i))

    def test_HESS_KnownGoodDipDers(self):
        self.longMessage = True
        self.assertEqual(self.oh.dipders.shape, self.dipders.shape)
        for i in range(self.oh.dipders.shape[0]):
            for j in range(self.oh.dipders.shape[1]):
                self.assertAlmostEqual(self.oh.dipders[i,j], \
                            self.dipders[i,j], delta=1e-6, \
                            msg="Dipole derivative element (" + str(i) + ',' + \
                                                        str(j) + ')')

    def test_HESS_KnownGoodIRSpectrum(self):
        self.longMessage = True
        self.assertEqual(self.oh.ir_comps.shape, self.ir_comps.shape)
        self.assertEqual(self.oh.ir_mags.shape, self.ir_mags.shape)
        for i in range(self.oh.ir_comps.shape[0]):
            self.assertAlmostEqual(self.oh.ir_mags[i,0], self.ir_mags[i,0], \
                    delta=1e-4, msg="IR T**2 element (" + str(i) + ')')
            for j in range(self.oh.ir_comps.shape[1]):
                self.assertAlmostEqual(self.oh.ir_comps[i,j], \
                            self.ir_comps[i,j], delta=1e-4, \
                            msg="IR T_i derivative element (" + str(i) + ',' + \
                                                        str(j) + ')')

    def test_HESS_KnownGoodPolDers(self):
        self.longMessage = True
        self.assertEqual(self.oh.polders.shape, self.polders.shape)
        for i in range(self.oh.polders.shape[0]):
            for j in range(self.oh.polders.shape[1]):
                self.assertAlmostEqual(self.oh.polders[i,j], \
                            self.polders[i,j], delta=1e-4, \
                            msg="Polarizability derivative element (" + \
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
                f.write(self.file_text_good \
                                    .replace(*self.bad_block_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any .hess files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in \
                                            self.bad_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_HESS_MissingBlockHess(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.hess_block, self.file_name + \
                    self.names.hess)

    def test_HESS_MissingBlockGeom(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.at_block, self.file_name + \
                    self.names.geom)

    def test_HESS_MissingBlockEnergy(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.energy, self.file_name + \
                    self.names.energy)

    def test_HESS_MissingBlockTemp(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.temp, self.file_name + \
                    self.names.temp)

    def test_HESS_MissingBlockFreqs(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.freq_block, self.file_name + \
                    self.names.freqs)

    def test_HESS_MissingBlockModes(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.modes_block, self.file_name + \
                    self.names.modes)

    def test_HESS_MissingBlockDipders(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name + \
                                            self.names.dipders).dipders)

    def test_HESS_MissingBlockIRSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name + \
                                            self.names.ir).ir_comps)
        self.assertIsNone(ORCA_HESS(self.file_name + \
                                            self.names.ir).ir_mags)

    def test_HESS_MissingBlockPolders(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        self.assertIsNone(ORCA_HESS(self.file_name + \
                                            self.names.polders).polders)


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
                f.write(self.file_text_good \
                                    .replace(*self.trunc_block_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + bname) for bname in \
                                            self.trunc_block_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_HESS_TruncatedBlocksHess(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.hess_block, self.file_name + self.names.hess)

    def test_HESS_TruncatedBlocksGeom(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.at_block, self.file_name + self.names.geom)

    def test_HESS_TruncatedBlocksFreqs(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.freq_block, self.file_name + self.names.freqs)

    def test_HESS_TruncatedBlocksModes(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.modes_block, self.file_name + self.names.modes)

    def test_HESS_TruncatedBlocksDipders(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.dipder_block, self.file_name + self.names.dipders)

    def test_HESS_TruncatedBlocksIRSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.ir_block, self.file_name + self.names.ir)

    def test_HESS_TruncatedBlocksPolders(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.polder_block, self.file_name + self.names.polders)

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
                f.write(self.file_text_good \
                                    .replace(*self.bad_data_substs[bname]))

    @classmethod
    def tearDownClass(self):
        # Remove any engrad files and try to remove the temp directory

        import os

        # Try to remove the files
        [os.remove(self.file_name + dname) for dname in \
                                            self.bad_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(self.testdir)

    def test_HESS_BadDataHessDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.hess_block, self.file_name + self.names.hess)

    def test_HESS_BadDataFreqDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.freq_block, self.file_name + self.names.freqs)

    def test_HESS_BadDataAtomSym(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        self.assertRaises(KeyError, ORCA_HESS, \
                                    self.file_name + self.names.atsym)

    def test_HESS_BadDataNormalModes(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # First dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.modes_block, self.file_name + self.names.modes)
        # Second dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.modes_block, \
                    self.file_name + self.names.modes + self.names.suffix_dim2)

    def test_HESS_BadDataDipdersDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.dipder_block, self.file_name + self.names.dipders)

    def test_HESS_BadDataIRSpectrum(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        # Bad dimension
        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.ir_block, self.file_name + self.names.ir)
        # Mismatched frequency
        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.ir_block, \
                    self.file_name + self.names.ir + self.names.suffix_badfreq)

    def test_HESS_BadDataPoldersDim(self):

        from opan.error import HESSError
        from opan.hess import ORCA_HESS

        assertErrorAndTypecode(self, HESSError, ORCA_HESS, \
                    HESSError.polder_block, self.file_name + self.names.polders)

## end class TestORCAHessBadData


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

    # Ensure correct typecode; suppress repeated error
    try:
        out = cobj(*args, **kwargs)
    except errtype as err:
        testclass.assertEqual(err.tc, tc)
    except:
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
    if os.path.isdir(dirname) or os.path.isfile(dirname):
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
