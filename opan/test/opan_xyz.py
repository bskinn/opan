#-------------------------------------------------------------------------------
# Name:        opan_xyz
# Purpose:     Test objects for opan.xyz
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     8 Mar 2016
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


class SuperOpanXYZ(unittest.TestCase):
    # Superclass for all OPAN .xyz test cases

    # Imports
    from textwrap import dedent
    import numpy as np
    from opan.test.utils import assertErrorAndTypecode

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

## end class SuperOpanXYZ


class TestOpanXYZGoodFileData(SuperOpanXYZ):
    # Ensuring importing a known OpenBabel xyz file with good data reports
    #  the correct geometric parameters, etc.

    @classmethod
    def setUpClass(cls):
        from opan.test.utils import setUpTestDir

        # Set up the directory and add the good file
        setUpTestDir(cls.testdir)

        # Write the file
        with open(cls.file_name, 'w') as f:
            f.write(cls.file_text_good)

    @classmethod
    def tearDownClass(cls):
        import os
        from opan.test.utils import tearDownTestDir

        # Delete the xyz file
        os.remove(cls.file_name)

        # Remove the test directory
        tearDownTestDir(cls.testdir)

    def setUp(self):
        # Load the object

        # Imports
        from opan.xyz import OpanXYZ

        # Create the object
        self.xyz = OpanXYZ(path=self.file_name)

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
        for t in zip(self.dist_Cu_O, self.xyz.dist_iter(None, 0, 1),
                                        range(self.dist_Cu_O.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-5,
                        msg="Cu-O distance mismatch at geom #" + str(t[2]) +
                        ": " + str(t[0:2]))
        for t in zip(self.dist_O_H1, self.xyz.dist_iter(None, 1, 2),
                                        range(self.dist_O_H1.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-5,
                        msg="O-H1 distance mismatch at geom #" + str(t[2]) +
                        ": " + str(t[0:2]))

    def test_XYZ_GoodFileDataAngles(self):
        for t in zip(self.angle_Cu_O_H1, self.xyz.angle_iter(None, 0, 1, 2),
                                        range(self.angle_Cu_O_H1.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-2,
                        msg="Cu-O-H1 angle mismatch at geom #" + str(t[2]) +
                        ": " + str(t[0:2]))

    def test_XYZ_GoodFileDataDihedrals(self):
        for t in zip(self.dihed_H2_O_Cu_H1,
                                    self.xyz.dihed_iter(None, 3, 1, 0, 2),
                                    range(self.dihed_H2_O_Cu_H1.shape[0])):
            self.assertAlmostEqual(t[0], t[1], delta=1e-2,
                        msg="H2-Cu-O-H1 dihedral mismatch at geom #" +
                        str(t[2]) + ": " + str(t[0:2]))

    def test_XYZ_GoodFileDataIterGeom(self):

        import numpy as np

        idxs = [1,4,8]
        for t in zip(np.array(self.geoms)[idxs],
                        self.xyz.geom_iter(idxs),
                        range(len(idxs))
                    ):
            for i in range(t[0].shape[0]):
                self.assertAlmostEqual(t[0][i], t[1][i], delta=1e-6,
                        msg="Geometry #" + str(t[2]) + \
                                ", coordinate element #" + str(i))

    #TEST: displ_iter call (dist_iter only calls displ_single)

## end class TestOpanXYZGoodData


class TestOpanXYZAltFileData(SuperOpanXYZ):
    # Ensuring successful import of an XYZ file with valid data of alternative
    #  formatting to that contained in SuperOpanXYZ.file_text_good
    #  (primarily if elements are specified by atomic number).

    @classmethod
    def setUpClass(cls):
        # Set up the directory and add alternate-g files
        from opan.test.utils import setUpTestDir
        setUpTestDir(cls.testdir)

        # Write the files
        for bname in cls.alt_file_data_substs.keys():
            with open(cls.file_name + bname, 'w') as f:
                f.write(cls.file_text_good
                                    .replace(*cls.alt_file_data_substs[bname]))

    @classmethod
    def tearDownClass(cls):
        # Remove any created files and try to remove the temp directory
        import os
        from opan.test.utils import tearDownTestDir

        # Try to remove the files
        [os.remove(cls.file_name + bname) for bname in
                                            cls.alt_file_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(cls.testdir)

    def test_XYZ_AltFileDataFirstCuAtomnum(self):

        from opan.xyz import OpanXYZ
        from opan.error import XYZError

        try:
            xyz = OpanXYZ(path=(self.file_name + self.names.Cu_1st_as_atomnum))
        except XYZError:  # pragma: no cover
            self.fail("XYZ import failed when success was expected.")

        self.assertEqual(xyz.atom_syms[0].upper(),
                                        self.good_direct_atoms[0])

    def test_XYZ_AltFileDataLateCuAtomnum(self):

        from opan.xyz import OpanXYZ
        from opan.error import XYZError

        try:
            xyz = OpanXYZ(path=(self.file_name +
                                            self.names.Cu_late_as_atomnum))
        except XYZError:  # pragma: no cover
            self.fail("XYZ import failed when success was expected.")

        self.assertEqual(xyz.atom_syms[0].upper(),
                                        self.good_direct_atoms[0])

## end class TestOpanXYZAltFileData


class TestOpanXYZGoodDirectData(SuperOpanXYZ):
    # Confirming sanity of an OpanXYZ generated directly from data.

    def setUp(self):
        # Load the object

        # Imports
        from opan.xyz import OpanXYZ

        # Create the object
        self.xyz = OpanXYZ(atom_syms=self.good_direct_atoms,
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
                    self.xyz.dist_iter(0,1,None),
                    range(self.good_direct_O1Dist.shape[0])
                        ):
            self.assertAlmostEqual(tup[0], tup[1], delta=1e-5,
                    msg="Distance between O1 and atom #" + str(tup[2]) +
                            " (" + self.xyz.atom_syms[tup[2]].capitalize() +
                            ")")

    def test_XYZ_GoodDirectDataIterCuO1Angle(self):
        for tup in zip(
                    self.good_direct_CuO1Angle,
                    self.xyz.angle_iter(0,0,1,None),
                    range(self.good_direct_CuO1Angle.shape[0])
                        ):
            self.assertAlmostEqual(tup[0], tup[1], delta=5e-3,
                    msg="Angle Cu-O-X with atom #" + str(tup[2]) +
                            " (" + self.xyz.atom_syms[tup[2]].capitalize() +
                            ")")

    def test_XYZ_GoodDirectDataIterDihed(self):
        for tup in zip(
                    self.good_direct_Dihed,
                    self.xyz.dihed_iter(0,None, 1,0,2),
                    range(self.good_direct_Dihed.shape[0])
                        ):
            self.assertAlmostEqual(tup[0], tup[1], delta=5e-3,
                    msg="Dihedral with atom #" + str(tup[2]) +
                            " (" + self.xyz.atom_syms[tup[2]].capitalize() +
                            ")")

## end class TestOpanXYZGoodDirectData


class TestOpanXYZBadFileData(SuperOpanXYZ):
    # Ensuring importing an XYZ file with data of generally valid formatting but
    #  invalid content raises the appropriate errors

    @classmethod
    def setUpClass(cls):
        # Set up the directory and add munged files
        from opan.test.utils import setUpTestDir
        setUpTestDir(cls.testdir)

        # Write the files
        for bname in cls.bad_file_data_substs.keys():
            with open(cls.file_name + bname, 'w') as f:
                f.write(cls.file_text_good
                                    .replace(*cls.bad_file_data_substs[bname]))

    @classmethod
    def tearDownClass(cls):
        # Remove any created files and try to remove the temp directory
        import os
        from opan.test.utils import tearDownTestDir

        # Try to remove the files
        [os.remove(cls.file_name + bname) for bname in
                                            cls.bad_file_data_substs.keys()]

        # Remove the directory
        tearDownTestDir(cls.testdir)

    def test_XYZ_BadFileDataFileStart(self):
        # File without a number-of-atoms spec at the very start should throw
        #  an XYZError

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.file_start))

    def test_XYZ_BadFileDataAtomCounts(self):
        # Multiple-geometry file without a consistent number of atoms
        #  throughout should throw an XYZError, regardless of the 'numats'
        #  spec at the top of each geometry block

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.atom_counts))

    def test_XYZ_BadFileDataConstNumats(self):
        # Multiple-geometry file should have identical num_ats specs at the
        #  head of each geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.const_numats))

    def test_XYZ_BadFileDataBadLateGeomAtomNum(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.bad_lategeom_atomnum))

    def test_XYZ_BadFileDataDifferentLateGeomAtomNum(self):
        # Discrepant atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                    XYZError.XYZFILE,
                    path=(self.file_name + self.names.diff_lategeom_atomnum))

    def test_XYZ_BadFileDataBadFirstGeomAtomNum(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.bad_1stgeom_atomnum))

    def test_XYZ_BadFileDataBadLateGeomAtomSym(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.bad_lategeom_atomsym))

    def test_XYZ_BadFileDataDifferentLateGeomAtomSym(self):
        # Discrepant atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                    XYZError.XYZFILE,
                    path=(self.file_name + self.names.diff_lategeom_atomsym))

    def test_XYZ_BadFileDataBadFirstGeomAtomSym(self):
        # Invalid atomic number in a non-initial geometry block.

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        self.assertErrorAndTypecode(XYZError, OpanXYZ,
                        XYZError.XYZFILE,
                        path=(self.file_name + self.names.bad_1stgeom_atomsym))

## end class TestOpanXYZBadFileData


class TestOpanXYZBadDirectData(SuperOpanXYZ):
    # Ensuring invalid data passed directly to construct an OpanXYZ instance
    #  throws the appropriate errors or otherwise exhibits the correct problems.

    def test_XYZ_BadDirectDataTruncCoords(self):
        from opan.xyz import OpanXYZ
        self.assertRaises(ValueError, OpanXYZ, \
                    coords=self.good_direct_geom[0][:-2], \
                    atom_syms=self.good_direct_atoms)

    def test_XYZ_BadDirectDataBadElement(self):
        from opan.xyz import OpanXYZ

        # Copy symbols, munge, and pass into assert
        munge_atoms = self.good_direct_atoms[:]
        munge_atoms[0] = 'CX'
        self.assertRaises(ValueError, OpanXYZ, \
                    coords=self.good_direct_geom, \
                    atom_syms=munge_atoms)

    def test_XYZ_BadDirectDataNonRealCoords(self):
        from opan.xyz import OpanXYZ
        import numpy as np

        # Copy coords, munge, and pass to assert
        munge_coords = self.good_direct_geom[0].copy() * 1.j
        self.assertRaises(ValueError, OpanXYZ, \
                    coords=munge_coords, \
                    atom_syms=self.good_direct_atoms)

    # TEST: ValueError for Angle at_1 == at_2 and at_3 == at_2 cases

    # TEST: ValueError for Dihed any at_x equal

    # TEST: IndexError for invalid (out-of-range) at_x for Dist, Angle, Dihed

    # TEST: XYZError.DIHED for too-nearly-linear atom trio(s)

    # TEST: ValueError when multiple 'None' values passed to an X_iter method

    # TEST: ValueError when 'None' passed to an X_iter with another iterable

## end class TestOpanXYZBadDirectData


class TestOpanXYZBadUsage(SuperOpanXYZ):
    # Ensuring bad use cases of OpanXYZ throw the appropriate errors or
    #  otherwise exhibit the correct problems.

    @classmethod
    def setUpClass(cls):
        # Set up the directory and add the good file
        from opan.test.utils import setUpTestDir
        setUpTestDir(cls.testdir)

        # Write the file
        with open(cls.file_name, 'w') as f:
            f.write(cls.file_text_good)

    @classmethod
    def tearDownClass(cls):
        import os
        from opan.test.utils import tearDownTestDir

        # Delete the xyz file
        os.remove(cls.file_name)

        # Remove the test directory
        tearDownTestDir(cls.testdir)

    def test_XYZ_BadUsageBadInitParams(self):
        # Calling without 'path' or both 'atom_syms' and 'coords'

        from opan.xyz import OpanXYZ

        self.assertRaises(NameError, OpanXYZ)
        self.assertRaises(NameError, OpanXYZ, atom_syms=['CU', 'O'])
        self.assertRaises(NameError, OpanXYZ, coords=[1,2,3])

    def test_XYZ_BadUsageReInitDirectData(self):
        # Should pitch a fit if someone uses the private data method to try to
        #  reinit an OpanXYZ instance

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        # Should work fine
        xyz = OpanXYZ(coords=self.good_direct_geom,
                                atom_syms=self.good_direct_atoms)

        # Should raise error to re-init with same contents.
        self.assertErrorAndTypecode(XYZError, xyz._load_data,
                        XYZError.OVERWRITE,
                        coords=self.good_direct_geom,
                        atom_syms=self.good_direct_atoms)

    def test_XYZ_BadUsageReInitFileData(self):
        # Should pitch a fit if someone uses the private file method to try to
        #  reinit an OpanXYZ instance

        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        # Should work fine
        xyz = OpanXYZ(path=self.file_name)

        # Should raise error to re-init with same contents.
        self.assertErrorAndTypecode(XYZError, xyz._load_file,
                        XYZError.OVERWRITE,
                        XYZ_path=self.file_name)

    def test_XYZ_BadUsageNotXYZFile(self):

        import os
        from opan.error import XYZError
        from opan.xyz import OpanXYZ

        badfilename = 'bad.file'

        try:
            with open(badfilename, 'w') as f:
                f.write("This is not an OpenBabel file.\n\n")
                f.write("(In case you were wondering...)\n\n")

            self.assertErrorAndTypecode(XYZError, OpanXYZ,
                            XYZError.XYZFILE,
                            path=badfilename)

        finally:
            os.remove(badfilename)


## end class TestOpanXYZBadUsage


# ========   TestSuite Generation Functions   ======= #

def suite_FileData():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanXYZAltFileData),
                tl.loadTestsFromTestCase(TestOpanXYZBadFileData),
                tl.loadTestsFromTestCase(TestOpanXYZGoodFileData)
                ])
    return s

## end def suite_FileData

def suite_DirectData():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanXYZBadDirectData),
                tl.loadTestsFromTestCase(TestOpanXYZGoodDirectData)
                ])
    return s

## end def suite_DirectData

def suite_Usage():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanXYZBadUsage)
                ])
    return s

## end def suite_Usage

def suite_All():
    s = unittest.TestSuite()
    s.addTests([suite_FileData(), suite_DirectData(), suite_Usage()])
    return s

## end def suite_All


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")
