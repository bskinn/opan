#-------------------------------------------------------------------------------
# Name:        utils.decorate
# Purpose:     Submodule containing custom decorators for OpenAnharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     30 Oct 2015
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

""" Custom decorators defined for OpenAnharmonic.

Decorators
----------

.. autoclass:: arraysqueeze

"""

# Imports
from functools import wraps as _wraps

# Decorators
class arraysqueeze(object):
    """ Converts selected arguments to squeezed np.arrays

    Pre-applies an ``np.array(...).squeeze()`` conversion to all positional
    arguments according to integer indices passed, and to any keyword arguments
    according to any strings passed.

    Each |int|_ argument passed instructs the decorator to convert the
    corresponding positional argument in the function definition.

    Each |str| argument passed instructs the decorator to convert the
    corresponding keyword argument.

    |str| parameters corresponding to keyword arguments absent in a
    particular function call and positional/optional argument indices
    beyond the range of the actual `\*args` of the decorated
    function are ignored.

    .. warning:: Likely fragile with optional arguments; needs to be tested.

    Arguments
    ---------
    \*args : |int|_ or |str|
        Arguments to convert to squeezed np.arrays.


    .. Decorator built using the class form per the exposition
        `here <http://www.artima.com/weblogs/viewpost.jsp?thread=240845>`__
        |extlink|.

    """

    def __init__(self, *args):
        """ Pass the positions of arguments to be arraysqueezed as integers
        """

        # Check for integers and strings
        for arg in args:
            if not (isinstance(arg, int) or isinstance(arg, str)):
                raise(ValueError("Invalid decorator argument: {0}".format(arg)))
            ## end if
        ## next arg

        # If all ok, store
        self.arglist = args

    ## end def __init__

    def __call__(self, f):
        """ Call the wrapped function after arraysqueezing selected arguments.

        Absent keyword arguments and positional/optional argument indices
        beyond the range of the actual *args of f are ignored.

        """
        @_wraps(f)
        def wrapped_f(*args, **kwargs):

            # Must import numpy
            import numpy as np

            # Working list of args values, since args is a tuple
            w_args = list(args)

            # Parse the arguments and arraysqueeze them. If squeezed to a
            #  singleton array, rewrap as a dimension-one array.
            for mod_arg in self.arglist:
                if isinstance(mod_arg, int) and mod_arg < len(w_args):
                    w_args[mod_arg] = np.array(w_args[mod_arg]).squeeze()
                    if not w_args[mod_arg].shape:
                        w_args[mod_arg] = w_args[mod_arg][np.newaxis]
                    ## end if
                elif isinstance(mod_arg, str) and mod_arg in kwargs:
                    kwargs[mod_arg] = np.array(kwargs[mod_arg]).squeeze()
                    if not kwargs[mod_arg].shape:
                        kwargs[mod_arg] = kwargs[mod_arg][np.newaxis]
                    ## end if
                # no 'else:' since type checked in __init__
                ## end if
            ## next mod_arg

            # Execute the function and return its result
            return f(*w_args, **kwargs)

        ## end def wrapped_f

        # Return the decorated function
        return wrapped_f

    ## end def __call__

## end class arraysqueeze




if __name__ == '__main__':
    print("Module not executable.")


