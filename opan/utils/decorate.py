#-------------------------------------------------------------------------------
# Name:        utils.decorate
# Purpose:     Submodule containing custom decorators for Open Anharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     30 Oct 2015
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

""" Custom decorators defined for Open Anharmonic.

Decorators
----------

.. autoclass:: arraysqueeze

.. autoclass:: kwarg_fetch

"""

# Imports
from functools import wraps as _wraps


# Decorators
class arraysqueeze(object):
    """ Converts selected arguments to squeezed np.arrays

    Pre-applies an ``np.array(...).squeeze()`` conversion to all positional
    arguments according to integer indices passed, and to any keyword arguments
    according to any strings passed.

    Each |int| argument passed instructs the decorator to convert the
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
    \*args : |int| or |str|
        Arguments to convert to squeezed |nparray|.


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
                raise ValueError("Invalid decorator argument: {0}".format(arg))

        # If all ok, store
        self.arglist = args

    def __call__(self, f):
        """ Call the wrapped function after arraysqueezing selected arguments.

        Absent keyword arguments and positional/optional argument indices
        beyond the range of the actual `*args` of `f` are ignored.

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
                    w_args[mod_arg] = np.asarray(w_args[mod_arg]).squeeze()
                    if not w_args[mod_arg].shape:
                        w_args[mod_arg] = w_args[mod_arg][np.newaxis]

                elif isinstance(mod_arg, str) and mod_arg in kwargs:
                    kwargs[mod_arg] = np.asarray(kwargs[mod_arg]).squeeze()
                    if not kwargs[mod_arg].shape:
                        kwargs[mod_arg] = kwargs[mod_arg][np.newaxis]

                # no 'else:' since type checked in __init__

            # Execute the function and return its result
            return f(*w_args, **kwargs)

        # end def wrapped_f

        # Return the decorated function
        return wrapped_f

    # end def __call__

# end class arraysqueeze


class kwarg_fetch(object):
    """ Fetch a function call to a missing or |None| keyword argument

    Arguments
    ---------
    (fill in here)



    .. warning::

        Fragile to needing to pass optional arguments to the callable?

    .. Decorator built using the class form per the exposition
        `here <http://www.artima.com/weblogs/viewpost.jsp?thread=240845>`__
        |extlink|.

    """

    import keyword as _kw

    @classmethod
    def ok_kwarg(cls, val):
        """ DOCSTRING """
        try:
            return str.isidentifier(val) and not cls._kw.iskeyword(val)
        except TypeError:
            # Catch integer values in particular; fail if found
            return False

    @classmethod
    def ok_tuplearg(cls, val):
        """ DOCSTRING """

        if not isinstance(val, tuple) and len(val) == 2 and \
                isinstance(val[0], int) and cls.ok_kwarg(val[1]):
            return False

    @classmethod
    def ok_argarg(cls, val):
        """ DOCSTRING """
        try:
            ok_kw = str.isidentifier(val) and not cls._kw.iskeyword(val)
        except TypeError:
            # Again, catch integer values in particular and fail if found
            ok_kw = False

        ok_pos = isinstance(val, int)

        ok_tup = cls.ok_tuplearg(val)

        return ok_kw or ok_pos or ok_tup

    def __init__(self, kw, c, *args):
        """Initialize with the keyword, callable, and relevant arguments

        """

        if not self.ok_kwarg(kw):
            raise ValueError("'kw' argument must be a valid non-keyword "
                             "identifier")
        self.kw = kw

        if not callable(c):
            raise ValueError("'c' argument must be callable")
        self.c = c

        if not all(map(self.ok_argarg, args)):
            raise ValueError("All 'args' must be valid non-keyword "
                             "identifier strings, integers, or (int, str) "
                             "tuples (valid-identifiers str's)")
        self.arglist = args

    def __call__(self, f):
        """Call the wrapped function after any needed fetch."""

        @_wraps(f)
        def wrapped_f(*args, **kwargs):
            # Working list of args, since args will be a tuple

            # Check for if the target kwarg is missing
            if self.kw not in kwargs or kwargs[self.kw] is None:
                # Missing. Must fetch.
                # Initialize as empty the arguments list to pass to the
                # callable
                arglist = []
                kwarglist = {}

                # Populate sequentially
                for arg in self.arglist:
                    # Argument is already parsed as ok, but this makes
                    # for a nice short way to check which type each argument
                    # is.
                    if self.ok_kwarg(arg):
                        # Keyword argument
                        kwarglist.update({arg: kwargs[arg]})
                    else:
                        # Positional argument
                        arglist.append(args[arg])

                # Call the callable and store the result into the target
                # keyword
                c_result = self.c(*arglist, **kwarglist)
                new_kwarg = {self.kw: c_result}
                kwargs.update(new_kwarg)

            # Either way, just call the wrapped function with the
            # originally-passed arguments, with the modified-or-not
            # kwargs
            return f(*args, **kwargs)

        # end def wrapped_f

        # Return the wrapped function
        return wrapped_f

    # end def __call__

# end class kwarg_fetch


if __name__ == '__main__':
    print("Module not executable.")


