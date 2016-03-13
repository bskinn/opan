#-------------------------------------------------------------------------------
# Name:        utils
# Purpose:     Helper functions for unit tests
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



def assertErrorAndTypecode(self, errtype, cobj, tc, *args, **kwargs):
    """ Wrapper for asserting correct OpanErrors and proper typecodes.

    Function tests (using testclass.assertX methods) whether 'cobj' raises
    'errtype' with typecode 'tc' when instantiated/called with *args and
    **kwargs.

    If imported into a :class:`~unittest.TestCase` subclass, this function
    is treated by Python as a class method and the initial `testclass`
    argument behaves like the `self` argument of such a class method.
    In particular, this means that the :class:`~unittest.TestCase` class
    should NOT be passed as the first argument in this usage situation.

    Parameters
    ----------
    self        : referenced object
        Subclass of unittest.TestCase (or related), from which the .assertX
        methods should be called
    errtype     : object reference
        Subclass of OpanError expected to be raised
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
    self.assertRaises(errtype, cobj, *args, **kwargs)

    # Ensure correct typecode; suppress repeated error, and ignore any
    #  other error raised, as it will have been reported by the above
    #  .assertRaises call
    try:
        out = cobj(*args, **kwargs)
    except errtype as err:
        self.assertEqual(err.tc, tc)
    except Exception:  # pragma: no cover
        pass

## end def assertErrorAndTypecode


def setUpTestDir(dirname):
    """ Create and change working directory to test directory.

    Parameters
    ----------
    dirname : str
        Name of desired working directory
    """

    import os, time

    # Wait 10ms for folder access to clear
    time.sleep(0.01)

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


if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")


