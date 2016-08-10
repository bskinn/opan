#-------------------------------------------------------------------------------
# Name:        opan_utils_decorate
# Purpose:     Test objects for opan.utils.decorate
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     12 Apr 2016
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


class TestOpanUtilsDecorateArraySqueeze(unittest.TestCase):

    from opan.utils.decorate import arraysqueeze as arsq

    @classmethod
    def setUpClass(cls):
        cls.longMessage = True

    @staticmethod
    @arsq(0)
    def fxn_2p_1d(p1, p2):
        return (p1, p2)

    @staticmethod
    @arsq(0, 1)
    def fxn_2p_2d(p1, p2):
        return (p1, p2)

    @staticmethod
    @arsq(1, 'test')
    def fxn_2p_1d_k(p1, p2, **kwargs):
        return p1, p2, kwargs.get('test')


    def test_arraysqueeze_GoodPartialWrap(self):
        import numpy as np

        # Store the result for detailed checking
        ret = self.fxn_2p_1d(1, 2)

        # First one needs to be an ndarray
        self.assertIsInstance(ret[0], np.ndarray,
                msg="First argument not converted to ndarray")
        self.assertTupleEqual(ret[0].shape, (1,),
                msg="Bad ndarray size for wrapped first argument")
        self.assertEqual(ret[0][0], 1,
                msg="Corrupted value in wrapped first argument")

        # Second one needs to still be an int
        self.assertIsInstance(ret[1], int,
                msg="Second argument not left as int")
        self.assertEqual(ret[1], 2,
                msg="Corrupted value in unmodified second argument")

    def test_arraysqueeze_GoodTwoArgWrap(self):
        import numpy as np

        # Store the result for detailed checking
        ret = self.fxn_2p_2d(1,2)

        # Both need to be ndarrays
        self.assertIsInstance(ret[0], np.ndarray,
                msg="First argument not converted to ndarray")
        self.assertTupleEqual(ret[0].shape, (1,),
                msg="Bad ndarray size for wrapped first argument")
        self.assertEqual(ret[0][0], 1,
                msg="Corrupted value in wrapped first argument")

        self.assertIsInstance(ret[1], np.ndarray,
                msg="Second argument not converted to ndarray")
        self.assertTupleEqual(ret[1].shape, (1,),
                msg="Bad ndarray size for wrapped second argument")
        self.assertEqual(ret[1][0], 2,
                msg="Corrupted value in wrapped second argument")

    def test_arraysqueeze_GoodArrayArgWrap(self):
        import numpy as np
        from opan.utils import pack_tups

        # Store the array result for detailed checking
        ary = np.array(range(1,9,2)).reshape((4,1))
        ret = self.fxn_2p_1d(ary,"x")

        # Confirm converted array matches expectations
        for i,t in enumerate(pack_tups(ary, ret[0])):
            self.assertEqual(t[0][0], t[1],
                msg="Mismatch at position {0} (Pre: {1}; Post: {2})"
                    .format(i, t[0][0], t[1]))

    def test_arraysqueeze_GoodKeywordArgWrap(self):
        import numpy as np

        # Result function with keyword
        ret = self.fxn_2p_1d_k(1, 2, test=3)

        # Confirm that the parameters are wrapped properly
        self.assertIsInstance(ret[1], np.ndarray,
                msg="Second argument not converted to ndarray")
        self.assertIsNotNone(ret[2],
                msg="Keyword argument improperly returned as None")
        self.assertIsInstance(ret[2], np.ndarray,
                msg="Keyword argument not converted to ndarray")

    def test_arraysqueeze_PosArgsPassedAsKwargs(self):

        import numpy as np

        ret = self.fxn_2p_1d_k(test=3, p2=2, p1=1)

        # Confirm that the parameters are wrapped properly
        self.assertIsInstance(ret[1], np.ndarray,
                              msg="Second argument not converted to ndarray")
        self.assertIsNotNone(ret[2],
                             msg="Keyword argument improperly returned as None")
        self.assertIsInstance(ret[2], np.ndarray,
                              msg="Keyword argument not converted to ndarray")

    def test_arraysqueeze_BadInitParam(self):

        with self.assertRaises(ValueError,
                    msg="Error not raised on float param index"):
            @self.arsq(1.5)
            def test_fxn(p1, p2):
                pass

        with self.assertRaises(ValueError,
                    msg="Error not raised on object param index"):
            @self.arsq(ValueError)
            def test_fxn(p1, p2):
                pass

        with self.assertRaises(ValueError,
                    msg="Error not raised on bool param index"):
            @self.arsq(True)
            def test_fxn(p1, p2):
                pass

        with self.assertRaises(ValueError,
                    msg="Error not raised on None param index"):
            @self.arsq(None)
            def test_fxn(p1, p2):
                pass

# end class TestOpanUtilsDecorateArraySqueeze


class TestOpanUtilsDecorateKwargFetch(unittest.TestCase):

    # These test functions are trivial, and only serve to probe
    #  fundamental behavior. Better intuitive sense can probably
    #  be made of kwarg_fetch by inspecting actual use cases in live code.

    from opan.utils.decorate import kwargfetch as kwf

    @classmethod
    def setUpClass(cls):
        cls.longMessage = True

    @staticmethod
    def f_1p(a):
        return a * a

    @staticmethod
    def f_2p(a, b):
        return b - a

    @staticmethod
    def f_1p_k(a, **kwargs):
        if 'kw' in kwargs and kwargs['kw']:
            return a + 2*kwargs['kw']
        else:
            return 2 * a

    @staticmethod
    def f_2p_1o(a, b, x=5):
        return a**b - a**x

    def test_kwargfetch_pos_to_pos(self):

        @self.kwf('c', self.f_1p, 0)
        def testfxn(a, **kwargs):
            return a * kwargs['c']

        self.assertEqual(testfxn(3), 27)
        # 'c' stored as f_1p(3) = 9.
        # Return value is then 3*9 = 27

    def test_kwargfetch_pos_to_opt_arg(self):

        @self.kwf('c', self.f_2p_1o, 0, 1, 3)
        def testfxn(a, b, x, y, **kwargs):
            return a - b - kwargs['c'] * (x + y)

        self.assertEqual(testfxn(3, 2, 3, 4), 505)
        # 'c' stored as f_2p_1o(3, 2, 4) = 9 - 81 = -72
        # Return value is then (3-2) - (-72)*(3+4) = 1 + 504 = 505

    def test_kwargfetch_pos_to_opt_kwarg(self):

        @self.kwf('c', self.f_2p_1o, 0, 1, x=3)
        def testfxn(a, b, x, y, **kwargs):
            return a - b - kwargs['c'] * (x + y)

        self.assertEqual(testfxn(3, 2, 3, 4), 505)
        # 'c' stored as f_2p_1o(3, 2, x=4) = 9 - 81 = -72
        # Return value is then (3-2) - (-72)*(7) = 1 + 504 = 505

    def test_kwargfetch_pos_to_kw(self):

        @self.kwf('c', self.f_1p_k, 1, kw=0)
        def testfxn(a, b, **kwargs):
            return (a + b)**kwargs['c']

        self.assertEqual(testfxn(5, 2), 13841287201)
        # 'c' stored as f_1p_k(2, kw=5) = 2 + 2*5 = 12
        # Return value is thus (5+2)**12 = 7**12 = 13841287201

    def test_kwargfetch_absent_opt_to_pos(self):

        @self.kwf('c', self.f_1p, 1)
        def testfxn(a, x=5, **kwargs):
            return (a + x) * (1 + kwargs['c'])

        self.assertEqual(testfxn(2), 182)
        # 'c' stored as f_1p(5) = 25
        # Return value is thus (2+5)*(1+25) = 7*26 = 182

    def test_kwargfetch_present_opt_to_pos(self):

        @self.kwf('c', self.f_1p, 1)
        def testfxn(a, x=5, **kwargs):
            return (a + x) * (1 + kwargs['c'])

        self.assertEqual(testfxn(2, 3), 50)
        # 'c' stored as f_1p(3) = 9
        # Return value is thus (2+3)*(1+9) = 5*10 = 50

    def test_kwargfetch_present_kwopt_to_pos(self):

        @self.kwf('c', self.f_1p, 1)
        def testfxn(a, x=5, **kwargs):
            return (a + x) * (1 + kwargs['c'])

        self.assertEqual(testfxn(2, x=3), 50)
        # 'c' stored as f_1p(3) = 9
        # Return value is thus (2+3)*(1+9) = 5*10 = 50

    def test_kwargfetch_absent_opt_to_opt_arg(self):

        @self.kwf('c', self.f_2p_1o, 0, 1, 2)
        def testfxn(a, b, x=-3, **kwargs):
            return (b + x) - (kwargs['c'] * a)

        self.assertEqual(testfxn(2, 8), -506.75)
        # 'c' stored as f_2p_1o(2, 8, -3) = 2**8 - 2**-3 = 256-0.125 = 255.875
        # Return value is thus (8-3) - (255.875*2) = 5 - 511.75 = -506.75

    def test_kwargfetch_absent_opt_to_opt_kwarg(self):

        @self.kwf('c', self.f_2p_1o, 1, 0, x=2)
        def testfxn(a, b, x=4, **kwargs):
            return (b + x) - (kwargs['c'] * a)

        self.assertEqual(testfxn(2, 3), 151)
        # 'c' stored as f_2p_1o(3, 2, x=4) = 3**2 - 3**4 = -72
        # Return value is thus (3-4) - (-72*2) = 7 + 144 = 151

    def test_kwargfetch_present_opt_to_opt_arg(self):

        @self.kwf('c', self.f_2p_1o, 1, 0, 2)
        def testfxn(a, b, x=13, **kwargs):
            return (b + x) - (kwargs['c'] * a)

        self.assertEqual(testfxn(1, 2, 3), 11)
        # 'c' stored as f_2p_1o(2, 1, 3) = 2**1 - 2**3 = 2-8 = -6
        # Return value is thus (2+3) - (1*-6) = 5+6 = 11

    def test_kwargfetch_present_kwopt_to_opt_arg(self):

        @self.kwf('c', self.f_2p_1o, 0, 1, 2)
        def testfxn(a, b, x=35, **kwargs):
            return (b + x) - (kwargs['c'] * a)

        self.assertEqual(testfxn(3, 5, x=4), -477)
        # 'c' stored as f_2p_1o(3, 5, x=4) = 3**5 - 3**4 = 243 - 81 = 162
        # Return value is thus (5+4) - (162*3) = 9 - 486 = -477

    def test_kwargfetch_present_opt_to_opt_kwarg(self):

        @self.kwf('c', self.f_2p_1o, 1, 0, x=2)
        def testfxn(a, b, x=15, **kwargs):
            return (b + x) - (kwargs['c'] * a)

        self.assertEqual(testfxn(5, 2, 3), -115)
        # 'c' stored as f_2p_1o(2, 5, x=3) = 2**5 - 2**3 = 32 - 8 = 24
        # Return value is (2+3) - (24*5) = 5 - 120 = -115

    def test_kwargfetch_present_kwopt_to_opt_kwarg(self):

        @self.kwf('c', self.f_2p_1o, 0, 1, x=2)
        def testfxn(a, b, x=-234234, **kwargs):
            return (b + x) - (kwargs['c'] * a)

        self.assertEqual(testfxn(3, 4, x=3), -155)
        # 'c' stored as f_2p_1o(3, 4, x=3) = 3**4 - 3**3 = 81 - 27 = 54
        # Return value is (4+3) - (54*3) = 7 - 162 = -155

    def test_kwargfetch_absent_opt_to_kw(self):

        @self.kwf('c', self.f_1p_k, 1, kw=2)
        def testfxn(a, b, x=4, **kwargs):
            return a * b * kwargs['c'] + x

        self.assertEqual(testfxn(5, 3), 169)
        # 'c' stored as f_1p_k(3, kw=4) = 3 + 2*4 = 11
        # Return value is 5*3*11 + 4 = 15*11 + 4 = 165 + 4 = 169

    def test_kwargfetch_present_opt_to_kw(self):

        @self.kwf('c', self.f_1p_k, 1, kw=2)
        def testfxn(a, b, x=5, **kwargs):
            return a * b + kwargs['c'] * x

        self.assertEqual(testfxn(2, 6, 3), 48)
        # 'c' stored as f_1p_k(6, kw=3) = 6 + 2*3 = 12
        # Return value is 2*6 + 12*3 = 12 + 36 = 48

    def test_kwargfetch_present_kwopt_to_kw(self):

        @self.kwf('c', self.f_1p_k, 0, kw=2)
        def testfxn(a, b, x=2, **kwargs):
            return a * b - kwargs['c'] * x

        self.assertEqual(testfxn(4, 2, x=3), -22)
        # 'c' stored as f_1p_k(4, kw=3) = 4 + 2*3 = 10
        # Return value is 4*2 - 10*3 = 8 - 30 = -22

    def test_kwargfetch_kw_to_pos(self):

        @self.kwf('c', self.f_1p, 'inv')
        def testfxn(**kwargs):
            return kwargs['c'] * (1 + kwargs['inv'])

        self.assertEqual(testfxn(inv=3), 36)
        # 'c' stored as f_1p(3) = 9
        # Return value is 9 * (1+3) = 36

    def test_kwargfetch_kw_to_opt_arg(self):

        @self.kwf('c', self.f_2p_1o, 0, 1, 'bip')
        def testfxn(a, b, **kwargs):
            return (a * kwargs['c']) + (b ** kwargs['bip'])

        self.assertEqual(testfxn(2, 4, bip=6), 4000)
        # 'c' stored as f_2p_1o(2, 4, 6) = 2**4 - 2**6 = 16 - 64 = -48
        # Return value is (2 * -48) + (4 ** 6) = -96 + 4096 = 4000

    def test_kwargfetch_kw_to_opt_kwarg(self):

        @self.kwf('c', self.f_2p_1o, 1, 0, x='dug')
        def testfxn(a, b, **kwargs):
            return (a + kwargs['c']) * (b + kwargs['dug'])

        self.assertEqual(testfxn(2, 5, dug=4), -5382)
        # 'c' stored as f_2p_1o(5,2,x=4) = 5**2 - 5**4 = 25-625 = -600
        # Return value is (2 + -600) * (5 + 4) = -598 * 9 = -5382

    def test_kwargfetch_kw_to_kw(self):

        @self.kwf('c', self.f_1p_k, 1, kw='bop')
        def testfxn(a, b, **kwargs):
            return (a + kwargs['c']) * (b + kwargs['bop'])

        self.assertEqual(testfxn(3, 1, bop=7), 144)
        # 'c' stored as f_1p_k(1, kw=7) = 1 + 14 = 15
        # Return value is (3+15) * (1+7) = 18*8 = 144

    def test_kwargfetch_absent_kw_to_kw(self):

        @self.kwf('c', self.f_1p_k, 1, kw='flip')
        def testfxn(a, b, **kwargs):
            return (a + kwargs['c']) * (b + kwargs['blip'])

        self.assertEqual(testfxn(2, 5, blip=7), 144)
        # 'c' stored as f_1p_k(5, kw=None) = 10
        # Return value is (2+10) * (5+7) = 12*12 = 144

    def test_kwargfetch_pos_to_kwpos(self):
        # Case where a positional arg is passed keyword-style
        @self.kwf('c', self.f_1p, 1)
        def testfxn(a, b, **kwargs):
            return (a + b) * kwargs['c']

        self.assertEqual(testfxn(3, b=5), 200)
        # 'c' stored as f_1p(5) = 25
        # Return value is (3+5) * 25 = 200

    def test_kwargfetch_allkwpos_to_kwpos(self):
        # Case where all positional args are passed keyword-style
        @self.kwf('c', self.f_1p, 1)
        def testfxn(a, b, **kwargs):
            return (a + b) * kwargs['c']

        self.assertEqual(testfxn(b=5, a=3), 200)
        # 'c' stored as f_1p(5) = 25
        # Return value is (3+5) * 25 = 200

    def test_kwargfetch_callablepos_as_kw(self):
        # Calling the callable with keyword arguments
        @self.kwf('c', self.f_1p, a=0)
        def testfxn(a, b, **kwargs):
            return (a + b) * kwargs['c']

        self.assertEqual(testfxn(3, 5), 72)
        # 'c' stored as f_1p(a=3) = 9
        # Return value is (3+5) * 9 = 72

    # Multiple kwarg_fetch decorators on a single function (simple, non-
    #  dependent case)
    def test_kwargfetch_multiple_independent(self):

        @self.kwf('c', self.f_1p, 0)
        @self.kwf('d', self.f_1p, 1)
        def testfxn(a, b, **kwargs):
            return (a + kwargs['d']) * (b + kwargs['c'])

        self.assertEqual(testfxn(3, 4), 247)
        # 'c' stored as f_1p(3) = 9
        # 'd' stored as f_1p(4) = 16
        # Return value is (3+16) * (4+9) = 19*13 = 247

    # Multiple kwarg_fetch decorators on a single function (complex case,
    #  where one kwarg_fetch depends on the kwarg injected by the other)
    def test_kwargfetch_multiple_dependent(self):

        @self.kwf('c', self.f_1p, 1)
        @self.kwf('d', self.f_1p, 'c')
        def testfxn(a, b, **kwargs):
            return (a + kwargs['d']) * (b + kwargs['c'])

        self.assertEqual(testfxn(5, 2), 126)
        # 'c' stored as f_1p(2) = 4
        # 'd' stored as f_1p(4) = 16
        # Return value is (5+16) * (2+4) = 21*6 = 126

    # Invalid target kwarg
    def test_kwargfetch_target_as_int(self):

        with self.assertRaises(ValueError, msg="Decorator init should "
                               "fail when integer passed for target keyword"):

            @self.kwf(1, self.f_1p, 0)
            def testfxn(a, **kwargs):
                pass

    def test_kwargfetch_nonidentifier_target(self):

        with self.assertRaises(ValueError, msg="Decorator init should "
                               "fail when target is not an identifier"):

            @self.kwf('abc@@#$', self.f_1p, 0)
            def testfxn(a, **kwargs):
                pass

    def test_kwargfetch_keyword_target(self):

        with self.assertRaises(ValueError, msg="Decorator init should "
                               "fail when target is a keyword"):

            @self.kwf('pass', self.f_1p, 0)
            def testfxn(a, **kwargs):
                pass

    # Invalid callable (not callable)
    def test_kwargfetch_noncallable_c(self):

        with self.assertRaises(TypeError, msg="Fetching callable "
                               "is not callable"):

            @self.kwf('c', 'x', 0)
            def testfxn(a, **kwargs):
                pass

    # Invalid arg list (e.g., too many or too few arg positions passed,
    #  invalid types passed to arg, indicated arg not present in
    #  wrapped function sig)
    def test_kwargfetch_too_few_posargs(self):

        @self.kwf('c', self.f_2p, 0)
        def testfxn(a, b, **kwargs):
            pass

        self.assertRaises(TypeError, testfxn, 1, 2,
                          msg="No TypeError raised when decorator called "
                              "with too few arguments for callable")

    def test_kwargfetch_too_many_posargs(self):

        @self.kwf('c', self.f_2p, 0, 1, 2)
        def testfxn(a, b, f, **kwargs):
            pass

        self.assertRaises(TypeError, testfxn, 1, 2, 3,
                          msg="No TypeError raised when decorator called "
                              "with too many arguments for callable")

    def test_kwargfetch_bad_posarg_index(self):

        @self.kwf('c', self.f_1p, 1)
        def testfxn(a, **kwargs):
            pass

        self.assertRaises(IndexError, testfxn, 1,
                          msg="No IndexError raised when decorator called "
                              "with positional argument not present in "
                              "wrapped function")

    # Invalid kwarg list
    # Python keyword as arg
    def test_kwargfetch_pykeyword_as_arg(self):

        with self.assertRaises(ValueError, msg="Decorator should fail when "
                                               "called with Python "
                                               "keyword as arg"):

            @self.kwf('c', self.f_1p_k, 'pass')
            def testfxn(a, **kwargs):
                pass

    # Python keyword as kwarg
    def test_kwargfetch_pykeyword_as_kwarg(self):

        with self.assertRaises(ValueError, msg="Decorator should fail when "
                                               "called with Python "
                                               "keyword as kwarg"):

            @self.kwf('c', self.f_1p_k, 0, bip='pass')
            def testfxn(a, **kwargs):
                pass

    # Invalid type passed to kw
    def test_kwargfetch_kwarg_invalid_type(self):

        with self.assertRaises(ValueError, msg="Decorator should fail when "
                                               "called with non-str/-int "
                                               "as kwarg"):
            @self.kwf('c', self.f_1p_k, 0, bip=IndexError())
            def testfxn(a, **kwargs):
                pass

    # Error if wrapped function doesn't allow kwargs
    def test_kwargfetch_no_kwargs_in_wrapped_function(self):

        @self.kwf('c', self.f_1p, 0)
        def testfxn(a):
            return 2*a

        self.assertRaises(TypeError, testfxn, 3,
                          msg="TypeError not raised when decorated function "
                              "without '**kwargs' is called")

    # Collision between injected kwarg and (optional-)positional name
    #  (if both positional and keyword specified for the same opt name)
    def test_kwargfetch_fetched_kwarg_collision_with_posarg(self):

        @self.kwf('a', self.f_1p, 0)
        def testfxn(a, **kwargs):
            pass

        @self.kwf('b', self.f_1p, 0)
        def testfxn2(a, b=3, **kwargs):
            pass

        self.assertRaises(TypeError, testfxn, 2,
                          msg="TypeError not raised when injected kwarg "
                              "collides with provided posarg")
        self.assertRaises(TypeError, testfxn2, 1, 2,
                          msg="TypeError not raised when injected kwarg "
                              "collides with provided posarg")

    # Decorator applied with no args at all
    def test_kwargfetch_no_args(self):

        with self.assertRaises(TypeError, msg="No TypeError raised "
                                              "when decorator applied "
                                              "with no arguments"):
            @self.kwf()
            def testfxn(a, **kwargs):
                pass

    # Decorator applied with only the target-keyword parameter
    def test_kwargfetch_only_kw_arg(self):

        with self.assertRaises(TypeError, msg="No TypeError raised when "
                                              "decorator applied with only "
                                              "a single argument"):
            @self.kwf('c')
            def testfxn(a, **kwargs):
                pass

    # Ensure no fetch performed if target kwarg already present
    def test_kwargfetch_no_fetch_if_kwarg_present(self):

        @self.kwf('c', self.f_1p, 0)
        def testfxn(a, **kwargs):
            return a + 2*kwargs['c']

        self.assertEquals(testfxn(3), 21)       # 'c' becomes 9 from f_1p
        self.assertEquals(testfxn(3, c=5), 13)  # 'c' not overridden

# end class TestOpanUtilsDecorateKwargFetch

def suite():
    s = unittest.TestSuite()
    tl = unittest.TestLoader()
    s.addTests([tl.loadTestsFromTestCase(TestOpanUtilsDecorateArraySqueeze),
                tl.loadTestsFromTestCase(TestOpanUtilsDecorateKwargFetch)
                ])
    return s


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

