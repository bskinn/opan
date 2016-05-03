#-------------------------------------------------------------------------------
# Name:        opan_utils_vector
# Purpose:     Test objects for opan.utils.vector
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     19 Apr 2016
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
from opan.test.utils import inject_tests


class TestOpanUtilsVectorParallelCheck(unittest.TestCase):

    def test_Utils_Vector_ParallelCheck_Good(self):
        from opan.utils.vector import parallel_check as pc
        import numpy as np

        # Parallel vectors
        self.assertTrue(pc(np.array([1, 2, 3]),
                           np.array([1.5, 3, 4.5])))

        # Anti-parallel vectors
        self.assertTrue(pc(np.array([-1, 5.3, 3.1]),
                           np.array([3, -15.9, -9.3])))

        # Non-(anti-)parallel vectors
        self.assertFalse(pc(np.array([4.8, 0.35, -1.822]),
                            np.array([-1.3, 3.77, 19.14])))

    def test_Utils_Vector_ParallelCheck_BadShape(self):
        from opan.utils.vector import parallel_check as pc
        import numpy as np

        self.assertRaises(ValueError, pc,
                          np.array([[1, 2, 3], [3, 2, 1]]),
                          np.array([2, 4, 9]))

    def test_Utils_Vector_ParallelCheck_LenMismatch(self):
        from opan.utils.vector import parallel_check as pc
        import numpy as np

        self.assertRaises(ValueError, pc,
                          np.array(range(3)),
                          np.array(range(4)))


class TestOpanUtilsVectorProjRejAngle(unittest.TestCase):
    import numpy as np
    from opan.const import OpanEnum

    # Enum for the types of data stored for calculations
    class DType(OpanEnum):
        V1 = 'V1'
        V2 = 'V2'
        PROJ = 'PROJ'
        REJ = 'REJ'
        ANG = 'ANG'

    # Enum for the types of vector pairs being tested
    class VecType(OpanEnum): # Types of vectors
        O1 = 'O1'   # Both order-one
        LOL = 'LOL' # Both large (large on large)
        SOS = 'SOS' # Both small (small on small)
        LOS = 'LOS' # Large onto small
        SOL = 'SOL' # Small onto large
        BS = 'BS'   # Badly-scaled

    # Enum for the relationship between the vectors being tested
    class RelType(OpanEnum): # Type of vector relationship
        NS = 'NS'      # Nonspecific
        PAR = 'PAR'    # Nearly parallel
        NORM = 'NORM'  # Nearly normal
        AP = 'AP'      # Nearly anti-parallel

    # Substitution string for naming the vector datasets
    namestr = "{0}_{1}"

    # Dict of dicts of data
    data = {
            # Unremarkable vectors with ~order-one components
            namestr.format(RelType.NS, VecType.O1):
                {DType.V1: np.array([1, 2, 3]),
                 DType.V2: np.array([-1, 3, 8]),
                 DType.PROJ: np.array([-0.391892, 1.175676, 3.135135]),
                 DType.REJ: np.array([1.391892, 0.824324, -0.135135]),
                 DType.ANG: np.float_(25.712002)},

            # Nearly-normal vectors with ~order-one components
            namestr.format(RelType.NORM, VecType.O1):
                {DType.V1: np.array([2, 8, -4, 2.5, -1.4]),
                 DType.V2: np.array([-1, 3, 6.2, 5, 7]),
                 DType.PROJ: np.array([0.000817, -0.002450, -0.005064,
                                       -0.004084, -0.005717]),
                 DType.REJ: np.array([1.999183, 8.002450, -3.994936,
                                      2.504084, -1.394283]),
                 DType.ANG: np.float_(90.053923)},

            # Nearly-parallel vectors with ~order-one components
            namestr.format(RelType.PAR, VecType.O1):
                {DType.V1: np.array([1, 2, 2.9999, 4]),
                 DType.V2: np.array([1.0001, 2, 3, 4]),
                 DType.PROJ: np.array([1.000087, 1.999973, 2.999960, 3.999947]),
                 DType.REJ: np.array([-0.000087, 0.000027,
                                      -0.000060, 0.000053]),
                 DType.ANG: np.float_(0.001267)},

            # Nearly-antiparallel vectors with ~order-one components
            namestr.format(RelType.AP, VecType.O1):
                {DType.V1: np.array([-0.5, 2.3, 1.4, -3.1]),
                 DType.V2: np.array([0.50001, -2.29999, -1.4, 3.1]),
                 DType.PROJ: np.array([-0.500011, 2.299992,
                                        1.400001, -3.100003]),
                 DType.REJ: np.array([0.000011, 0.000008,
                                     -0.000001, 0.000003]),
                 DType.ANG: np.float_(179.999814)},

            # Two large vectors far from parallel/normal
            namestr.format(RelType.NS, VecType.LOL):
                {DType.V1: np.array([376328., 384874.,
                                     992834., 182873.]),
                 DType.V2: np.array([538344., 283747.,
                                     1838447., 929292.]),
                 DType.PROJ: np.array([269185.658799, 141880.699195,
                                       919270.144855, 464669.577884]),
                 DType.REJ: np.array([107142.341201, 242993.300805,
                                      73563.855145, -281796.577884]),
                 DType.ANG: np.float_(20.151554)},

            # Two small vectors far from parallel/normal
            namestr.format(RelType.NS, VecType.SOS):
                {DType.V1: np.array([0.000045, -0.00031,
                                     -0.000915, 0.000002]),
                 DType.V2: np.array([0.0002874, -0.0003987,
                                     0.0000034, 0.000719]),
                 DType.PROJ: np.array([0.000051, -0.000071,
                                       0.000001, 0.000128]),
                 DType.REJ: np.array([-0.000006, -0.000239,
                                      -0.000916, -0.000126]),
                 DType.ANG: np.float_(80.787151)},

            # Large onto small, far from parallel/normal
            namestr.format(RelType.NS, VecType.LOS):
                {DType.V1: np.array([238973., 239884.,
                                     -1092938., -893983.]),
                 DType.V2: np.array([0.0002874, -0.0003987,
                                     0.0000034, 0.000719]),
                 DType.PROJ: np.array([-255163.218951, 353979.037564,
                                       -3018.632375, -638351.963904]),
                 DType.REJ: np.array([494136.218951, -114095.037564,
                                      -1089919.367625, -255631.036096]),
                 DType.ANG: np.float_(122.176632)},

            # Small onto large, far from parallel/normal
            namestr.format(RelType.NS, VecType.SOL):
                {DType.V1: np.array([0.00000374, -0.0000233,
                                     0.0002837, 0.0000026]),
                 DType.V2: np.array([538344., 283747., 1838447., 929292.]),
                 DType.PROJ: np.array([0.000061, 0.000032,
                                       0.000207, 0.000105]),
                 DType.REJ: np.array([-0.000057, -0.000055,
                                       0.000077, -0.000102]),
                 DType.ANG: np.float_(31.859107)},

            # Two badly scaled vectors far from parallel/normal
            namestr.format(RelType.NS, VecType.BS):
                {DType.V1: np.array([0.000015, 6214., -0.000235, 12374.]),
                 DType.V2: np.array([-0.00005, 38184., 0.000045, 21669.]),
                 DType.PROJ: np.array([-0.000013, 10011.853795,
                                       0.000012, 5681.616904]),
                 DType.REJ: np.array([0.000028, -3797.853795,
                                      -0.000247, 6692.383096]),
                 DType.ANG: np.float_(33.760587)},


            # Two large vectors nearly parallel
            namestr.format(RelType.PAR, VecType.LOL):
                {DType.V1: np.array([554387., 38185., -532247., 12374.]),
                 DType.V2: np.array([554389., 38184., -532248., 12375.]),
                 DType.PROJ: np.array([554387.488030, 38183.895862,
                                    -532246.548415, 12374.966250]),
                 DType.REJ: np.array([-0.488030, 1.104138,
                                      -0.451585, -0.966250]),
                 DType.ANG: np.float_(0.000120)},

            # Two small vectors nearly parallel
            namestr.format(RelType.PAR, VecType.SOS):
                {DType.V1: np.array([0.000015, 0.000016, -0.000042,
                                    -0.000034, 0.000065, -0.000033]),
                 DType.V2: np.array([0.000014, 0.000017, -0.000041,
                                    -0.000033, 0.000066, -0.000032]),
                 DType.PROJ: np.array([0.000014, 0.000017, -0.000041,
                                      -0.000033, 0.000066, -0.000032]),
                 DType.REJ: np.array([0.000001, -0.000001, -0.000001,
                                     -0.000001, -0.000001, -0.000001]),
                 DType.ANG: np.float_(1.483536)},

            # Large onto small, nearly parallel
            namestr.format(RelType.PAR, VecType.LOS):
                {DType.V1: np.array([14001., 17002., -41003., -33001.,
                                     66004., -32005., 75008.]),
                 DType.V2: np.array([0.000014, 0.000017, -0.000041,
                                    -0.000033, 0.000066, -0.000032, 0.000075]),
                 DType.PROJ: np.array([14001.205610, 17001.463955,
                                      -41003.530715, -33002.841795,
                                       66005.683590, -32002.755680,
                                       75006.458626]),
                 DType.REJ: np.array([-0.205610, 0.536045, 0.530715,
                                       1.841795, -1.683590, -2.244320,
                                       1.541374]),
                 DType.ANG: np.float_(0.001811)},

            # Small onto large, nearly parallel
            namestr.format(RelType.PAR, VecType.SOL):
                {DType.V1: np.array([1.20E-05, -1.30E-05, 3.80E-05,
                                    -7.90E-05, -8.50E-05, 2.20E-05,
                                     4.50E-05]),
                 DType.V2: np.array([12006., -13009., 38007.,
                                    -79008., -85006., 22001.,
                                     45003.]),
                 DType.PROJ: np.array([0.000012, -0.000013, 0.000038,
                                      -0.000079, -0.000085, 0.000022,
                                       0.000045]),
                 DType.REJ: np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                 DType.ANG: np.float_(0.004356)},

            # Two badly scaled vectors nearly parallel
            namestr.format(RelType.PAR, VecType.BS):
                {DType.V1: np.array([0.00053, -2135., 0.00015, 65548.,
                                    -0.0085, 0.00022, -17815., -0.00035]),
                 DType.V2: np.array([0.00054, -2136., 0.00014, 65549.,
                                    -0.0084, 0.00023, -17816., -0.00036]),
                 DType.PROJ: np.array([0.000540, -2135.960458, 0.000140,
                                       65547.786547, -0.008400, 0.000230,
                                      -17815.670188, -0.000360]),
                 DType.REJ: np.array([-0.000010, 0.960458, 0.000010,
                                       0.213453, -0.000100, -0.000010,
                                       0.670188,0.000010]),
                 DType.ANG: np.float_(0.001004)},


            # Two large vectors nearly normal
            namestr.format(RelType.NORM, VecType.LOL):
                {DType.V1: np.array([654564., 48249., -248896.,
                                    -54789., 24444.]),
                 DType.V2: np.array([-6048., 68589., -34061.,
                                      95549., -108100.]),
                 DType.PROJ: np.array([11.146008, -126.404358, 62.771856,
                                     -176.089606, 199.220153]),
                 DType.REJ: np.array([654552.853992, 48375.404358,
                                     -248958.771856, -54612.910394,
                                       24244.779847]),
                 DType.ANG: np.float_(90.024498)},

            # Two small vectors nearly normal
            namestr.format(RelType.NORM, VecType.SOS):
                {DType.V1: np.array([0.000048, 0.000017, -0.000032,
                                    0.000091, -0.000016]),
                 DType.V2: np.array([0.000015, 0.000007, 0.00001,
                                    0.000001, 0.000035]),
                 DType.PROJ: np.array([0.000000, 0.000000, 0.000000,
                                    0.000000, 0.000001]),
                 DType.REJ: np.array([0.000048, 0.000017, -0.000032,
                                    0.000091, -0.000017]),
                 DType.ANG: np.float_(89.350346)},

            # Large onto small, nearly normal
            namestr.format(RelType.NORM, VecType.LOS):
                {DType.V1: np.array([32185., 27265., -30185.,
                                    108115., -13755.]),
                 DType.V2: np.array([0.000015, 0.000007, 0.00001,
                                    0.000001, 0.000035]),
                 DType.PROJ: np.array([-14.343750, -6.693750, -9.562500,
                                        -0.956250, -33.468750]),
                 DType.REJ: np.array([32199.343750, 27271.693750,
                                    -30175.437500, 108115.956250,
                                    -13721.531250]),
                 DType.ANG: np.float_(90.018157)},

            # Small onto large, nearly normal
            namestr.format(RelType.NORM, VecType.SOL):
                {DType.V1: np.array([0.000028, -0.000035, 0.000017,
                                     0.000098, 0.000072, -0.000055]),
                DType.V2: np.array([16589., 22852., 44185.,
                                    14273., 24599., 65489.]),
                 DType.PROJ: np.array([0., 0., 0., 0., 0., 0.]),
                 DType.REJ: np.array([0.000028, -0.000035, 0.000017,
                                      0.000098, 0.000072, -0.000055]),
                 DType.ANG: np.float_(90.073867)},

            # Two badly scaled vectors nearly normal
            namestr.format(RelType.NORM, VecType.BS):
                {DType.V1: np.array([0.0015, 6214., 2319., 145.]),
                 DType.V2: np.array([8285., -0.0004, 0.0034, 2166.]),
                 DType.PROJ: np.array([35.485053, -0.000002,
                                       0.000015, 9.277082]),
                 DType.REJ: np.array([-35.483553, 6214.000002,
                                      2318.999985, 135.722918]),
                 DType.ANG: np.float_(89.683234)},


            # Two large vectors nearly anti-parallel
            namestr.format(RelType.AP, VecType.LOL):
                {DType.V1: np.array([215484., 665452., -654587.,
                                     541887., 64657., -6546347.,
                                    -687887., -1137889.]),
                 DType.V2: np.array([-215485., -665453., 654589.,
                                     -541882., -64659., 6546378.,
                                      687889., 1137888.]),
                 DType.PROJ: np.array([215484.046715, 665450.056100,
                                      -654586.104161, 541879.602766,
                                       64658.713955, -6546349.039449,
                                      -687885.956845, -1137882.966092]),
                 DType.REJ: np.array([-0.046715, 1.943900, -0.895839,
                                       7.397234, -1.713955, 2.039449,
                                      -1.043155, -6.033908]),
                 DType.ANG: np.float_(179.999914)},

            # Two small vectors nearly anti-parallel
            namestr.format(RelType.AP, VecType.SOS):
                {DType.V1: np.array([0.000041, -0.000038, 0.000091,
                                    -0.000019, 0.000037, -0.000068,
                                    -0.000071, -0.000055]),
                 DType.V2: np.array([-0.000041, 0.000039, -0.00009,
                                      0.000019, -0.000036, 0.000068,
                                      0.00007, 0.000055]),
                 DType.PROJ: np.array([0.000041, -0.000039, 0.000091,
                                      -0.000019, 0.000036, -0.000068,
                                      -0.000070, -0.000055]),
                 DType.REJ: np.array([0.000000, 0.000001, 0.000000,
                                      0.000000, 0.000001, 0.000000,
                                     -0.000001, 0.000000]),
                 DType.ANG: np.float_(179.379006)},

            # Large onto small, nearly anti-parallel
            namestr.format(RelType.AP, VecType.LOS):
                {DType.V1: np.array([41002., -38997., 90004., -18997.,
                                     36001., -68002., -70003., -54988.]),
                 DType.V2: np.array([-0.000041, 0.000039, -0.00009,
                                      0.000019, -0.000036, 0.000068,
                                      0.00007, 0.000055]),
                 DType.PROJ: np.array([40999.983927, -38999.984711,
                                       89999.964717, -18999.992551,
                                       35999.985887, -67999.973342,
                                      -69999.972558, -54999.978438]),
                 DType.REJ: np.array([2.016073, 2.984711, 4.035283,
                                      2.992551, 1.014113, -2.026658,
                                     -3.027442, 11.978438]),
                 DType.ANG: np.float_(179.994978)},

            # Small onto large, nearly anti-parallel
            namestr.format(RelType.AP, VecType.SOL):
                {DType.V1: np.array([0.000038, -0.000044, 0.000057,
                                     0.000098, 0.000089, -0.000072,
                                    -0.000022, 0.000025, -0.000017]),
                 DType.V2: np.array([-38007., 44009., -56995.,
                                     -98004., -88999., 71995.,
                                      22002., -25001., 16996.]),
                 DType.PROJ: np.array([0.000038, -0.000044, 0.000057,
                                       0.000098, 0.000089, -0.000072,
                                      -0.000022, 0.000025, -0.000017]),
                 DType.REJ: np.array([0., 0., 0., 0., 0., 0., 0., 0., 0.]),
                 DType.ANG: np.float_(179.995212)},

            # Two badly scaled vectors nearly anti-parallel
            namestr.format(RelType.AP, VecType.BS):
                {DType.V1: np.array([25778., -35778., 0.000032, -47789.,
                                     -0.000038, 0.000041, 24448., -35779.,
                                     -0.000017]),
                 DType.V2: np.array([-25779., 35772., -0.000031, 47788.,
                                      0.000038, -0.000041, -24444., 35777.,
                                      0.000016]),
                 DType.PROJ: np.array([25780.714146, -35774.378619, 0.000031,
                                      -47791.177610, -0.000038, 0.000041,
                                       24445.625376, -35779.378952, -0.000016]),
                 DType.REJ: np.array([-2.714146, -3.621381, 0.000001,
                                       2.177610, 0.000000, 0.000000,
                                       2.374624, 0.378952, -0.000001]),
                 DType.ANG: np.float_(179.995917)}
            }

    # Template functions
    # Vector projection template
    def template_proj(self, name, data):
        from opan.utils.vector import proj

        v1 = data[self.DType.V1]
        v2 = data[self.DType.V2]
        p = proj(v1, v2)
        for i, t in enumerate(zip(p, data[self.DType.PROJ])):
            self.assertAlmostEqual(*t, delta=1e-6,
                       msg="Test {0}: Index {1}; V1 = {2}; V2 = {3}"
                            .format(name, i, v1, v2))

    # Vector rejection template
    def template_rej(self, name, data):
        from opan.utils.vector import rej

        v1 = data[self.DType.V1]
        v2 = data[self.DType.V2]
        r = rej(v1, v2)
        for i, t in enumerate(zip(r, data[self.DType.REJ])):
            self.assertAlmostEqual(*t, delta=1e-6,
                       msg="Test {0}: Index {1}; V1 = {2}; V2 = {3}"
                            .format(name, i, v1, v2))

    # Vector angle template
    def template_angle(self, name, data):
        from opan.utils.vector import vec_angle

        v1 = data[self.DType.V1]
        v2 = data[self.DType.V2]
        a = vec_angle(v1, v2)
        self.assertAlmostEqual(a, data[self.DType.ANG], delta=1e-6,
                   msg="Test {0}: V1 = {1}; V2 = {2}"
                   .format(name, v1, v2))

    # Populate the local namespace with the auto-generated
    #  test methods
    _locals = locals()
    inject_tests(_locals, data, "test_Vector_Proj_Good_{0}", template_proj)
    inject_tests(_locals, data, "test_Vector_Rej_Good_{0}", template_rej)
    inject_tests(_locals, data, "test_Vector_Angle_Good_{0}", template_angle)


    def setUp(self):
        self.longMessage = True


    def test_Utils_Vector_Proj_BadVec_NotVector(self):
        import numpy as np
        from opan.utils.vector import proj

        self.assertRaises(ValueError, proj,
                          np.array(range(16)).reshape((4, 4)),
                          np.array(range(16)))

    def test_Utils_Vector_Proj_BadVecOnto_NotVector(self):
        import numpy as np
        from opan.utils.vector import proj

        self.assertRaises(ValueError, proj,
                          np.array(range(16)),
                          np.array(range(16)).reshape((4, 4)))

    def test_Utils_Vector_Proj_BadVecsShapeMismatch(self):
        import numpy as np
        from opan.utils.vector import proj

        self.assertRaises(ValueError, proj,
                          np.array(range(5)), np.array(range(6)))


class TestOpanUtilsVectorOrthoBasis(unittest.TestCase):

    def setUp(self):
        self.longMessage = True


def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanUtilsVectorParallelCheck),
                tl.loadTestsFromTestCase(TestOpanUtilsVectorProjRejAngle),
                tl.loadTestsFromTestCase(TestOpanUtilsVectorOrthoBasis)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")
