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
    class blocknames(object):
        numats = 'numats'
        energy = 'energy'
        grad = 'grad'
        geom = 'geom'
        E = frozenset([numats, energy, grad, geom])

    bad_block_substs = {
            blocknames.numats: ('umber of', 'asdlkjf'),
            blocknames.energy: ('urrent total en', 'alksdjfdsflkj'),
            blocknames.grad: ('Eh/bohr', 'asdlkjsdf'),
            blocknames.geom: ('inates in Bohr', 'asldkjfas;ldkfj')
                        }

    trunc_block_substs = {
            blocknames.grad: ('-0.000007912155', '#-0.000007912155'),
            blocknames.geom: ('1     2.6798241', '#1     2.6798241')
                            }

    bad_atom_subst = (' 29    -1.5545432', '229    -1.5545432')

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

## end class SuperORCAEngrad


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
                        self.coords[i,0], delta=1e-7, \
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

        import os

        # Check if test directory already exists (or file of same name);
        #  error if so
        if os.path.isdir(self.testdir) or os.path.isfile(self.testdir):
            raise(IOError("Cannot create new test directory!"))

        # Create and change to test directory
        os.mkdir(self.testdir)
        os.chdir(self.testdir)

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

        # Switch to parent directory
        os.chdir(os.path.pardir)

        # Try to remove the temp directory
        os.rmdir(self.testdir)

    def test_ENGRAD_MissingBlockNumAts(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.numats, self.file_name + self.blocknames.numats)

    def test_ENGRAD_MissingBlockEnergy(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.en, self.file_name + self.blocknames.energy)

    def test_ENGRAD_MissingBlockGrad(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.gradblock, self.file_name + self.blocknames.grad)

    def test_ENGRAD_MissingBlockGeom(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.geomblock, self.file_name + self.blocknames.geom)

## end def TestORCAEngradMissingBlocks


class TestORCAEngradTruncatedBlocks(SuperORCAEngrad):
    # Ensuring importing a file with incomplete grad or geom block throws
    #  the right errors

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

        # Switch to parent directory
        os.chdir(os.path.pardir)

        # Try to remove the temp directory
        os.rmdir(self.testdir)

    def test_ENGRAD_TruncatedBlockGrad(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.gradblock, self.file_name + self.blocknames.grad)

    def test_ENGRAD_TruncatedBlockGeom(self):

        from opan.error import GRADError
        from opan.grad import ORCA_ENGRAD

        assertErrorAndTypecode(self, GRADError, ORCA_ENGRAD, \
                    GRADError.geomblock, self.file_name + self.blocknames.geom)

## end class TestORCAEngradTruncatedBlocks

#TODO: Test for invalid atomic number in geom block

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
      0

    $act_coord
      0

    $act_energy
            0.000000

    $hessian
    12
                      0          1          2          3          4          5
          0       0.066469   0.015972   0.015971  -0.064711  -0.013251  -0.013250
          1       0.015973   0.012022   0.001164  -0.008930  -0.019230   0.002723
          2       0.015971   0.001164   0.012021  -0.008929   0.002723  -0.019230
          3      -0.064442  -0.008981  -0.008980   0.279514   0.126830   0.126832
          4      -0.013261  -0.019222   0.002957   0.126794   0.386720  -0.186574
          5      -0.013260   0.002957  -0.019221   0.126795  -0.186574   0.386722
          6      -0.000879  -0.008875   0.001833  -0.107561  -0.144984   0.031446
          7      -0.001710   0.003676  -0.001257  -0.132983  -0.302391   0.105565
          8      -0.001010  -0.002629   0.003531   0.015123   0.078038  -0.065118
          9      -0.000879   0.001832  -0.008876  -0.107566   0.031445  -0.144988
         10      -0.001011   0.003531  -0.002629   0.015122  -0.065115   0.078033
         11      -0.001710  -0.001257   0.003677  -0.132987   0.105559  -0.302390
                      6          7          8          9         10         11
          0      -0.000879  -0.001710  -0.001010  -0.000879  -0.001011  -0.001710
          1      -0.008875   0.003677  -0.002630   0.001833   0.003531  -0.001257
          2       0.001833  -0.001257   0.003531  -0.008875  -0.002630   0.003677
          3      -0.107534  -0.132994   0.015146  -0.107539   0.015144  -0.132998
          4      -0.144996  -0.302394   0.078052   0.031463  -0.065104   0.105565
          5       0.031465   0.105571  -0.065107  -0.145000   0.078046  -0.302394
          6       0.098093   0.139953  -0.028019   0.010347   0.013907  -0.005260
          7       0.139953   0.304719  -0.106612  -0.005259  -0.006004   0.002304
          8      -0.028019  -0.106612   0.067592   0.013906   0.031202  -0.006005
          9       0.010347  -0.005259   0.013906   0.098098  -0.028018   0.139957
         10       0.013907  -0.006004   0.031202  -0.028018   0.067589  -0.106606
         11      -0.005260   0.002304  -0.006005   0.139957  -0.106606   0.304718

    $vibrational_frequencies
    12
        0        0.000000
        1        0.000000
        2        0.000000
        3        0.000000
        4        0.000000
        5        0.000000
        6      374.347011
        7      654.278526
        8      741.107742
        9     1588.036185
       10     3348.605257
       11     3401.267490

    $normal_modes
    12 12
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
                      6          7          8          9         10         11
          0      -0.147481   0.003956  -0.000000   0.000449  -0.000028  -0.000000
          1      -0.038968  -0.004174   0.007852  -0.000264  -0.000005   0.000074
          2      -0.038964  -0.004174  -0.007852  -0.000264  -0.000005  -0.000074
          3       0.515899  -0.072606  -0.000001  -0.057431   0.036190  -0.000000
          4       0.140132   0.064410  -0.058388  -0.030030   0.024191  -0.051126
          5       0.140116   0.064405   0.058395  -0.030029   0.024192   0.051126
          6       0.554841   0.451533  -0.622852   0.441612  -0.286292   0.298539
          7       0.114660  -0.355151   0.128248  -0.054015  -0.606060   0.606872
          8       0.117915  -0.403979  -0.303516   0.547272   0.222412  -0.199971
          9       0.554835   0.451475   0.622882   0.441609  -0.286316  -0.298536
         10       0.117926  -0.404000   0.303485   0.547279   0.222406   0.199950
         11       0.114648  -0.355136  -0.128290  -0.054029  -0.606085  -0.606845

    #
    # The atoms: label  mass x y z
    #
    $atoms
    4
     Cu    63.5500     -1.554543    -0.492302    -0.492219
     O     15.9990      1.864044     0.366037     0.366022
     H      1.0080      2.679824     1.989221    -0.162252
     H      1.0080      2.679852    -0.162202     1.989204

    $actual_temperature
      0.000000

    $dipole_derivatives
    12
        -0.423570    -0.050796    -0.050789
        -0.086092    -0.148223    -0.007875
        -0.086082    -0.007875    -0.148219
        -0.731090    -0.035803    -0.035804
        -0.034041    -0.810897     0.095446
        -0.034047     0.095445    -0.810898
         0.577331     0.103902    -0.017308
         0.289558     0.582050     0.002537
        -0.169432    -0.090117     0.377080
         0.577336    -0.017302     0.103901
        -0.169426     0.377074    -0.090108
         0.289561     0.002548     0.582041

    #
    # The IR spectrum
    #  wavenumber T**2 TX TY  TY
    #
    $ir_spectrum
    12
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
          0.00       0.0000       0.0000       0.0000       0.0000
        374.35      17.2149       4.0855       0.5115       0.5114
        654.28     373.0665      14.1743      -9.2780      -9.2776
        741.11      64.5960       0.0004       5.6827      -5.6836
       1588.04     165.1818      10.1060       5.6148       5.6146
       3348.61     782.7819     -23.9756     -10.1969     -10.1968
       3401.27     499.6110       0.0003      15.8054     -15.8051


    $end
    """)
    testdir = 'orca_test_dir'
    file_name = 'test.hess'

    # Matching matrix values
    hess = np.matrix(
       [[ 0.066469,  0.015972,  0.015971, -0.064711, -0.013251, -0.01325 ,
        -0.000879, -0.00171 , -0.00101 , -0.000879, -0.001011, -0.00171 ],
       [ 0.015973,  0.012022,  0.001164, -0.00893 , -0.01923 ,  0.002723,
        -0.008875,  0.003677, -0.00263 ,  0.001833,  0.003531, -0.001257],
       [ 0.015971,  0.001164,  0.012021, -0.008929,  0.002723, -0.01923 ,
         0.001833, -0.001257,  0.003531, -0.008875, -0.00263 ,  0.003677],
       [-0.064442, -0.008981, -0.00898 ,  0.279514,  0.12683 ,  0.126832,
        -0.107534, -0.132994,  0.015146, -0.107539,  0.015144, -0.132998],
       [-0.013261, -0.019222,  0.002957,  0.126794,  0.38672 , -0.186574,
        -0.144996, -0.302394,  0.078052,  0.031463, -0.065104,  0.105565],
       [-0.01326 ,  0.002957, -0.019221,  0.126795, -0.186574,  0.386722,
         0.031465,  0.105571, -0.065107, -0.145   ,  0.078046, -0.302394],
       [-0.000879, -0.008875,  0.001833, -0.107561, -0.144984,  0.031446,
         0.098093,  0.139953, -0.028019,  0.010347,  0.013907, -0.00526 ],
       [-0.00171 ,  0.003676, -0.001257, -0.132983, -0.302391,  0.105565,
         0.139953,  0.304719, -0.106612, -0.005259, -0.006004,  0.002304],
       [-0.00101 , -0.002629,  0.003531,  0.015123,  0.078038, -0.065118,
        -0.028019, -0.106612,  0.067592,  0.013906,  0.031202, -0.006005],
       [-0.000879,  0.001832, -0.008876, -0.107566,  0.031445, -0.144988,
         0.010347, -0.005259,  0.013906,  0.098098, -0.028018,  0.139957],
       [-0.001011,  0.003531, -0.002629,  0.015122, -0.065115,  0.078033,
         0.013907, -0.006004,  0.031202, -0.028018,  0.067589, -0.106606],
       [-0.00171 , -0.001257,  0.003677, -0.132987,  0.105559, -0.30239 ,
        -0.00526 ,  0.002304, -0.006005,  0.139957, -0.106606,  0.304718]])

    atoms = np.array([['CU'],['O'],['H'],['H']])
    coords = np.matrix([-1.5545432, -0.4923021, -0.4922193, \
                            1.8640436, 0.3660366, 0.3660223, \
                            2.6798241, 1.9892213, -0.1622520, \
                            2.6798521, -0.1622023, 1.9892043]).transpose()
    masses = np.matrix([63.55, 15.999, 1.008, 1.008]).transpose()

## end class SuperORCAHess


class TestORCAHessKnownGood(SuperORCAHess):
    # Testing errors related to file parsing of a known-good HESS file

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

        # Delete the hess file
        os.remove(self.file_name)

        # Switch to parent directory
        os.chdir(os.path.pardir)

        # Try to remove the temp directory
        os.rmdir(self.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.hess import ORCA_HESS

        # Create the object
        self.oh = ORCA_HESS(self.file_name)

    def test_HESS_KnownGoodAtomVec(self):
        # Confirm the values coming out of the HESS match the known-good
        #  example file.

        # Confirm the atom list is good
        for i in range(self.oh.atom_syms.shape[0]):
            self.assertEqual(self.oh.atom_syms[i,0], self.atoms[i,0])

    def test_HESS_KnownGoodAtomMasses(self):
        # Confirm the values coming out of the HESS match the known-good
        #  example file.

        # Confirm the atom masses are good
        self.longMessage = True
        for i in range(self.oh.atom_masses.shape[0]):
            self.assertAlmostEqual(self.oh.atom_masses[i,0], \
                        self.masses[i,0], delta=1e-4, \
                        msg="Atom index " + str(i))

    def test_HESS_KnownGoodHess(self):
        # Confirm the known-good Hessian matches what's expected.

        # Confirm the Hessian matrix is good
        self.longMessage = True
        for i in range(self.oh.hess.shape[0]):
            for j in range(self.oh.hess.shape[1]):
                self.assertAlmostEqual(self.oh.hess[i,j], \
                            self.hess[i,j], delta=1e-6, \
                            msg="Hessian element (" + str(i) + ',' + \
                                                        str(j) + ')')

    def test_HESS_KnownGoodGeom(self):
        self.longMessage = True
        for i in range(self.oh.geom.shape[0]):
            self.assertAlmostEqual(self.oh.geom[i,0], \
                        self.coords[i,0], delta=1e-6, \
                        msg="Coordinate index " + str(i))

    def test_HESS_KnownGoodInitFlagDefined(self):
        self.assertTrue('initialized' in self.oh.__dict__)
        if 'initialized' in self.oh.__dict__:
            self.assertTrue(self.oh.initialized)

## end class TestORCAEngradKnownGood

#RESUME: MissingBlock, TruncBlock, invalid atom tests

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


if __name__ == '__main__':

    # Run test package, default verbose. Can override with '-q' at
    #  commandline
    unittest.main(verbosity=2)

## end main block
